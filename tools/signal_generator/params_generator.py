import numpy as np
import os


def write_params_txt(output_path: str,
                     fft_eq: np.ndarray,
                     blank_table: np.ndarray,
                     blank_range: int,
                     blank_valid_size: int,
                     mf_coeff: np.ndarray,
                     power_shift: int,
                     comp_mode: int,
                     comp_width: int,
                     scale_comp: int,
                     iq_reverse: int) -> None:

    with open(output_path, 'w') as f:
        # fft_eq: 4096 个 uint16，逗号分隔
        f.write("fft_eq=")
        f.write(','.join(str(x) for x in fft_eq))
        f.write('\n')

        # blank_table: 31 个索引值
        f.write("blank_table=")
        f.write(','.join(str(x) for x in blank_table))
        f.write('\n')

        # 标量参数
        f.write(f"blank_range={blank_range}\n")
        f.write(f"blank_valid_size={blank_valid_size}\n")

        # mf_coeff: 8 个系数
        f.write("mf_coeff=")
        f.write(','.join(str(x) for x in mf_coeff))
        f.write('\n')

        f.write(f"power_shift={power_shift}\n")
        f.write(f"comp_mode={comp_mode}\n")
        f.write(f"comp_width={comp_width}\n")
        f.write(f"scale_comp={scale_comp}\n")
        f.write(f"iq_reverse={iq_reverse}\n")

    print(f"Generated: {output_path}")


def default() -> None:
    """默认参数配置"""
    fft_eq = np.ones(4096, np.uint16) * 0x8000
    blank_table = np.array([2036,2037,2038,2039,2040,2041,2042,2043,2044,2045,2046,2047,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    blank_range=4
    blank_valid_size=5
    mf_coeff=np.array([7085,6144,4007,1965,725,0,0,0])
    power_shift=5
    comp_mode=2
    comp_width=4
    scale_comp=6964
    iq_reverse=0

    # 生成到默认位置
    output_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'app', 'IQ_mismatch_data')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'params.txt')

    write_params_txt(output_path, fft_eq, blank_table, blank_range, blank_valid_size,
                     mf_coeff, power_shift, comp_mode, comp_width, scale_comp, iq_reverse)


def main() -> None:

    fft_eq = np.ones(4096, np.uint16) * 0x8000
    blank_table = np.array([2036,2037,2038,2039,2040,2041,2042,2043,2044,2045,2046,2047,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    blank_range=4
    blank_valid_size=5
    mf_coeff=np.array([7085,6144,4007,1965,725,0,0,0])
    power_shift=5
    comp_mode=2
    comp_width=4
    scale_comp=6964
    iq_reverse=0

    output_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'app', 'IQ_mismatch_data')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'params.txt')

    write_params_txt(output_path, fft_eq, blank_table, blank_range, blank_valid_size,
                     mf_coeff, power_shift, comp_mode, comp_width, scale_comp, iq_reverse)


if __name__ == "__main__":
    main()