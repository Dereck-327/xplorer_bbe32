
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class IQSignal :
    def __init__(self, leak_Q_I: float = 0.15, # Q对I的串扰
                 leak_I_Q: float = 0.10, # I对Q的串扰
                 gain_imblance: float = 1.05, # IQ 增益不平衡
                 noise_floor: float = -20,
                 fs: int = 1e9, f_sig: float = 50e6 #% 信号频率默认50MHz
                 , N: int = 4096, seed: int = 42) :
        self.N = N
        Fs = fs  # 采样率 默认1GHz
        Ts = 1/Fs
        t = np.linspace(0, N, N) * Ts
        sig_I = np.cos(2 * np.pi * f_sig * t)
        sig_Q = np.sin(2 * np.pi * f_sig * t)
        self._f_axis = np.linspace(-N/2, N/2, N) * Fs/N

        # 模拟硬件混合
        I_raw = gain_imblance * (sig_I + leak_Q_I * sig_Q)
        Q_raw = 1.0 * (sig_Q + leak_I_Q * sig_I)
        # 底噪
        rng = np.random.default_rng(seed=seed)  # 可选设置 seed 以保证结果可复现
        self._I_raw = I_raw + 10 ** (noise_floor / 20) * rng.random(I_raw.size)
        self._Q_raw = Q_raw + 10 ** (noise_floor / 20) * rng.random(Q_raw.size)
        # 去除直流
        self._I_ac = self._I_raw - np.mean(self._I_raw)
        self._Q_ac = self._Q_raw - np.mean(self._Q_raw)

        self._rho = self.calculate_rho()
        self._gain_corr = self.calculate_g_corr()
        pass

    def calculate_rho(self) -> float:
        """ 计算相关系数 rho
            公式: I_new = I, Q_new = Q - rho * I
                使得 E[I_new * Q_new] = 0

        """
        E_II = np.mean(self._I_ac ** 2)
        E_IQ = np.mean(self._I_ac * self._Q_ac)
        rho = E_IQ / E_II
        print(f"相关系数rho :: {rho}")
        return rho

    def calculate_g_corr(self) -> float:
        """ 计算增益平衡系数 gain_corr
            使得 E[I^2] = E[Q^2]

        """
        I_orth = self._I_ac
        Q_orth = self._Q_ac - self._rho * self._I_ac
        pwr_I = np.mean(I_orth ** 2)
        pwr_Q = np.mean(Q_orth ** 2)
        gain_corr = np.sqrt(pwr_I / pwr_Q)
        print(f"增益平衡系数 :: {gain_corr}")
        return gain_corr

    def gram_schmidt(self) -> None:
        """ 正交化的IQ失配盲补偿算法

        """
        I_F = np.fft.fft(self._I_ac)
        Q_F = np.fft.fft(self._Q_ac)
        Q_F_corr = Q_F - self._rho * I_F  # 去相关
        Q_F_corr = self._gain_corr * Q_F_corr  # 增益平衡
        I_F_corr = I_F

        sig_cal_F = I_F_corr + 1j * Q_F_corr
        sig_cal_F = np.fft.fftshift(np.abs(sig_cal_F) ** 2) / self.N  # 取幅度谱

        sig_cal_F = 10 * np.log10(sig_cal_F)
        self.sig_cal_F = sig_cal_F - np.max(sig_cal_F)
        pass

    def get_distorted_spectrum(self) -> np.array:
        sig_distroted = self._I_raw + 1j * self._Q_raw
        pwr_raw = np.fft.fftshift(np.abs(np.fft.fft(sig_distroted)) ** 2) / self.N
        pwr_raw_db = 10 * np.log10(pwr_raw)
        pwr_raw_db = pwr_raw_db - np.max(pwr_raw_db)
        return pwr_raw_db

    def generate_fft_iq_ac(self, path: str) -> None:
        """ 将 I_F=fft(I_ac)、Q_F=fft(Q_ac) 定点化 (S16 块浮点) 写文件。
            复数 -> 交织存 real/imag; 全部分量共用一个块指数 fft_bexp。
            真值 = mantissa * 2^(fft_bexp - 15)
        """
        I_F = np.fft.fft(self._I_ac)
        Q_F = np.fft.fft(self._Q_ac)

        # 交织成 [re, im, re, im, ...]
        i_iq = np.empty(2 * self.N); i_iq[0::2] = I_F.real; i_iq[1::2] = I_F.imag
        q_iq = np.empty(2 * self.N); q_iq[0::2] = Q_F.real; q_iq[1::2] = Q_F.imag

        # 两路共用一个块指数, 让最大幅值缩进 S16 满量程 (<1.0)
        peak = float(np.max(np.abs(np.concatenate([i_iq, q_iq]))))
        fft_bexp = int(np.ceil(np.log2(peak))) if peak > 0 else 0
        scale = 2.0 ** fft_bexp

        def quant(x):
            m = np.round(x / scale * (1 << 15)).astype(np.int64)
            return np.clip(m, -(1 << 15), (1 << 15) - 1).astype(np.int16)

        i_mant = quant(i_iq)
        q_mant = quant(q_iq)

        with open(path, "w") as f:
            f.write("i_f=" + ",".join(str(int(v)) for v in i_mant) + "\n")
            f.write("q_f=" + ",".join(str(int(v)) for v in q_mant) + "\n")
            f.write(f"fft_bexp={fft_bexp}\n")


def main() -> None:
    import os

    out_dir = os.path.join(os.path.dirname(__file__), "out")
    os.makedirs(out_dir, exist_ok=True)

    # 生成 CHIRP_NUM=4 组数据 (不同噪声种子)
    for n in range(4):
        seed=42 + n
        iq_sig = IQSignal(seed=seed)

        path = os.path.join(out_dir, f"ra_data_{n}.txt")
        iq_sig.generate_fft_iq_ac(path)
        # print(f"写出 {path}")

    # 可视化最后一组的校准前后频谱
    iq_sig.gram_schmidt()
    f_axis = iq_sig._f_axis / 1e6  # MHz
    plt.subplot(2, 1, 1)
    plt.plot(f_axis, iq_sig.get_distorted_spectrum(), 'r', linewidth=1.5)
    plt.title('校准前频谱 (存在镜像串扰)')
    plt.ylabel('归一化功率 (dB)')
    plt.ylim([-100, 10]); plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(f_axis, iq_sig.sig_cal_F, 'b', linewidth=1.5)
    plt.title('校准后频谱 (镜像被抑制)')
    plt.xlabel('频率 (MHz)'); plt.ylabel('归一化功率 (dB)')
    plt.ylim([-100, 10]); plt.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__" :
    main()