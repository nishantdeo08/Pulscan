import glob
import numpy as np
import re

def parse_par_file(par_file_path):
    """Parses F0 (spin freq), A1 (proj. semi-major axis), M2 (companion mass) from par file."""
    f0 = a1 = m2 = None
    with open(par_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == 'F0':
                f0 = float(parts[1])
            elif parts[0] == 'A1':
                a1 = float(parts[1])
            elif parts[0] == 'M2':
                m2 = float(parts[1])
    if f0 is None:
        raise ValueError("F0 not found in .par file")
    return f0, a1, m2

def matching_frequency(par_file_path, pulscan_output_path):
    DM_array = []
    sigma_array = []

    matching_files = glob.glob(pulscan_output_path)
    time = 10 * 60  # seconds

    # Get F0 (spin frequency), A1, M2 from par file
    freq_expect, A1, M2 = parse_par_file(par_file_path)

    for file in matching_files:
        data = np.loadtxt(file, skiprows=1, delimiter=',')
        if data.ndim == 1:
            data = data.reshape(1, -1)  # If only one candidate in file

        r = data[:, 2]
        sigma = data[:, 0]
        freq_achieve = r / time

        for i, freq in enumerate(freq_achieve):
            if 0.99 * freq_expect <= freq <= 1.01 * freq_expect:
                match = re.search(r'DM(\d+\.\d+)', file)
                if match:
                    DM_val = float(match.group(1))
                    DM_array.append(DM_val)
                    sigma_array.append(sigma[i])

    return DM_array, sigma_array

par_file_path = '/lustre_archive/spotlight/Nishant/Par_files/Simulation1010.par'
pulscan_output_path = '/lustre_archive/spotlight/Nishant/data/Simulation1010_DM*.gpucand'

DMs, sigmas = matching_frequency(par_file_path, pulscan_output_path)






