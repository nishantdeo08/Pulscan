import os
import glob
import numpy as np
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def extract_par_parameters(par_file_path):
    """
    Extracts key parameters from a .par file:
    - F0: Spin frequency
    - A1: Projected semi-major axis
    - M2: Companion mass
    - DM: Dispersion Measure
    Returns these values as floats.
    """
    f0 = a1 = m2 = dm = None
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
            elif parts[0] == 'DM':
                dm = float(parts[1])

    if f0 is None:
        raise ValueError("F0 not found in .par file")
    return f0, a1, m2, dm

def extract_spectrum_time(spectrum_file_path):
    """
    Extracts:
    - Number of bins in time series
    - Width of each time series bin (seconds)
    from the .spectrum file.
    """
    num_bins = None
    samp_time = None

    with open(spectrum_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Number of bins in the time series'):
                num_bins = float(line.split('=')[-1].strip())
            elif line.startswith('Width of each time series bin (sec)'):
                samp_time = float(line.split('=')[-1].strip())

    if num_bins is None or samp_time is None:
        raise ValueError(f"Missing values in spectrum file: {spectrum_file_path}")

    return num_bins, samp_time


def extract_matching_candidates(par_file_path, pulscan_output_path, spectrum_file_path):
    """
    Matches expected F0 from .par file with candidate frequencies in .gpucand files.
    Returns arrays of matching DM values and their corresponding SNRs (sigmas).
    """
    DM_array = []
    sigma_array = []

    matching_files = glob.glob(pulscan_output_path)
    num_bins, samp_time = extract_spectrum_time(spectrum_file_path)
    time = num_bins * samp_time  # Integration time in seconds

    freq_expect, A1, M2, DM = extract_par_parameters(par_file_path)

    for file in matching_files:
        try:
            data = np.genfromtxt(file, delimiter=',', skip_header=1)

            # Skip empty or invalid files
            if data.size == 0:
                raise ValueError("Empty or invalid data")

            if data.ndim == 1:
                data = data.reshape(1, -1)

            if data.shape[1] < 3:
                raise ValueError("Not enough columns")

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

        except Exception as e:
            print(f"Skipping {file} due to error: {e}")
            continue

        #print("DM Array: ", DM_array)
        #print("Sigma Array: ", sigma_array)

    return np.array(DM_array), np.array(sigma_array)


def check_pulsar_candidate(par_file_path, DMs, sigmas):
    """
    Checks if the DM from the .par file is found in the matched DM list.
    Returns (DM, sigma) if match is found, else returns False.
    """
    _, _, _, DM = extract_par_parameters(par_file_path)
    DM_round = round(DM, 1)

    np.set_printoptions(threshold=np.inf)

    print(DM_round)

    print(np.sort(DMs))

    if DM_round in DMs:
        idx = np.where(DMs == DM_round)[0]
        print("Correct")
        return DM_round, sigmas[idx]
    else:
        print("Wrong")
        return False


def plot_candidate_result(par_file_path, DMs, sigmas, ax):
    """
    Adds a 3D scatter point to the shared figure.
    - X = A1
    - Y = M2
    - Z = sigma (SNR)
    """
    _, A1, M2, _ = extract_par_parameters(par_file_path)
    result = check_pulsar_candidate(par_file_path, DMs, sigmas)

    if result is not False:
        _, sigma = result
        ax.scatter(A1, M2, sigma, color='blue', marker='o', label='Pulsar Candidate')
    else:
        ax.scatter(A1, M2, 0, color='red', marker='x', label='No Match')

    print(f"Done for {os.path.basename(par_file_path)}")


def main():
    """
    Main driver:
    - Loops over all .par files
    - Matches corresponding .gpucand candidate files
    - Parses parameters, finds matches, and visualizes results
    """
    par_file = '/lustre_archive/spotlight/Nishant/Par_files/Simulation1010.par'
    # par_files = glob.glob(par_file_pattern)

    spectrum_file_path = '/lustre_archive/spotlight/Nishant/dedispersed/Simulation1010_DM10.00.inf'

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    base_name = os.path.basename(par_file)
    sim_id = os.path.splitext(base_name)[0]
    pulscan_output_path = f'/lustre_archive/spotlight/Nishant/data/{sim_id}_DM*.gpucand'

    DMs, sigmas = extract_matching_candidates(par_file, pulscan_output_path, spectrum_file_path)
    plot_candidate_result(par_file, DMs, sigmas, ax)

    # Set plot labels
    ax.set_xlabel('Projected Semi-Major Axis (A1)')
    ax.set_ylabel('Companion Mass (M2)')
    ax.set_zlabel('SNR')
    ax.set_title("Pulsar Candidate Summary")

    # Optional: de-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    '''
    # Save the figure
    output_dir = "./plots"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pulsar_candidates_3Dplot.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure to: {output_path}")
    '''
    # Show the plot
    plt.show()

if __name__ == "__main__":
    main()
