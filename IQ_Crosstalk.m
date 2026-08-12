clear; clc; close all;

%% 1. 系统参数设置
Fs = 1e9;               % 采样率 1 GHz
Ts = 1/Fs;
f_sig = 50e6;           % 信号频率 50 MHz
N = 4096;              % 采样点数 (2^14)，保证频谱分辨率
t = (0:N-1) * Ts;

% 模拟 ADC 硬件缺陷参数 (隔离度差)
% 理想情况：I与Q完全隔离
% 实际情况：I路信号会漏到Q路，Q路也会漏到I路
leak_Q_to_I = 0.15;     % Q 对 I 的串扰 (比如 -16dB 隔离度)
leak_I_to_Q = 0.10;     % I 对 Q 的串扰 (不对称串扰)
gain_imbalance = 1.05;  % 增益不平衡 (I路比Q路大一点)

%% 2. 信号生成
% 生成理想的单音复信号: exp(j*2*pi*f*t)
% 理想情况下，频谱只应在 +50MHz 处有峰值
sig_ideal_I = cos(2*pi*f_sig*t)';
sig_ideal_Q = sin(2*pi*f_sig*t)';

% 模拟硬件混合 (物理层的串扰)
% I_adc = Gain_I * (I_ideal + k1 * Q_ideal)
% Q_adc = Gain_Q * (Q_ideal + k2 * I_ideal)
I_raw = gain_imbalance * (sig_ideal_I + leak_Q_to_I * sig_ideal_Q);
Q_raw = 1.0 * (sig_ideal_Q + leak_I_to_Q * sig_ideal_I);

% 添加底噪 (模拟真实 ADC SNR)
noise_floor = -20; % dB
I_raw = I_raw + 10^(noise_floor/20) * randn(size(I_raw));
Q_raw = Q_raw + 10^(noise_floor/20) * randn(size(Q_raw));

% 组合成未校准的复信号
sig_distorted = I_raw + 1j*Q_raw;

%% 3. IQ 盲校准算法 (Gram-Schmidt 正交化)
% ---------------------------------------------------------
% 步骤 A: 去除直流 (DC Offset)
I_ac = I_raw - mean(I_raw);
Q_ac = Q_raw - mean(Q_raw);

% 步骤 B: 相位/串扰校正 (通过去相关)
% 计算相关系数 rho
% 公式: I_new = I, Q_new = Q - rho * I
% 使得 E[I_new * Q_new] = 0
E_II = mean(I_ac.^2);
E_IQ = mean(I_ac .* Q_ac);
rho = E_IQ / E_II;

I_orth = I_ac;
Q_orth = Q_ac - rho * I_ac;

% 步骤 C: 增益平衡 (Gain Correction)
% 使得 E[I^2] = E[Q^2]
pwr_I = mean(I_orth.^2);
pwr_Q = mean(Q_orth.^2);
g_corr = sqrt(pwr_I / pwr_Q);

I_cal = I_orth;
Q_cal = Q_orth * g_corr;

% 组合成校准后的复信号
sig_calibrated = I_cal + 1j*Q_cal;

%% 4. 功率谱密度 (PSD) 分析与可视化
nfft = 4096;
window = hamming(nfft);
overlap = nfft/2;

Pxx_raw = fftshift(abs(fft(sig_distorted)).^2)/N; % 归一化功率谱
Pxx_cal = fftshift(abs(fft(sig_calibrated)).^2)/N; % 归一化功率谱
f_axis = (-N/2:N/2-1)*(Fs/N); % 频率轴

% 转换为 dBFS (归一化到主峰)
Pxx_raw_dB = 10*log10(Pxx_raw);
Pxx_raw_dB = Pxx_raw_dB - max(Pxx_raw_dB); % 归一化

Pxx_cal_dB = 10*log10(Pxx_cal);
Pxx_cal_dB = Pxx_cal_dB - max(Pxx_cal_dB); % 归一化

% 绘图
figure('Position', [100, 100, 1000, 600]);

subplot(3,1,1);
plot(f_axis/1e6, Pxx_raw_dB, 'r', 'LineWidth', 1.5);
grid on;
title(['校准前频谱 (存在严重的镜像串扰) - 信号: 50MHz, 采样: 1GHz']);
xlabel('频率 (MHz)'); ylabel('归一化功率 (dB)');
xlim([-200, 200]); ylim([-100, 10]);
% 标注镜像
text(-50, -10, '\leftarrow 镜像干扰 (Image)', 'Color', 'r', 'FontSize', 12);
text(50, 5, '主信号 \rightarrow', 'HorizontalAlignment', 'right', 'FontSize', 12);

subplot(3,1,2);
plot(f_axis/1e6, Pxx_cal_dB, 'b', 'LineWidth', 1.5);
grid on;
title('时域校准后频谱 (镜像被抑制)');
xlabel('频率 (MHz)'); ylabel('归一化功率 (dB)');
xlim([-200, 200]); ylim([-100, 10]);
text(-50, -80, '\leftarrow 镜像被压制到底噪附近', 'Color', 'b', 'FontSize', 12);

%% 频域处理计算
I_F = fft(I_ac);
Q_F = fft(Q_ac);

Q_F_corr = Q_F - rho * I_F;   % 去相关
Q_F_corr = g_corr * Q_F_corr; % 增益平衡

I_F_corr = I_F;
sig_cal_F = I_F_corr + 1j*Q_F_corr;
sig_cal_F = fftshift(abs(sig_cal_F).^2)/N; % 取幅度谱

sig_cal_F = 10*log10(sig_cal_F);
sig_cal_F = sig_cal_F - max(sig_cal_F); % 归一化

subplot(3,1,3);
plot(f_axis/1e6, sig_cal_F, 'b', 'LineWidth', 1.5);
grid on;
title('频域校准后频谱 (镜像被抑制)');
xlabel('频率 (MHz)'); ylabel('归一化功率 (dB)');
xlim([-200, 200]); ylim([-100, 10]);
text(-50, -80, '\leftarrow 镜像被压制到底噪附近', 'Color', 'b', 'FontSize', 12);

