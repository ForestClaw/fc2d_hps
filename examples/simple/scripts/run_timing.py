import sys
import os
import subprocess
import numpy as np

class SimpleResults:
    def __init__(self):
        self.L1_error = 0
        self.L2_error = 0
        self.LI_error = 0
        self.setup_time = 0
        self.build_time = 0
        self.solve_time = 0
        self.copy_time = 0
        self.n_levels = 0
        self.n_cells = 0
        self.dofs = 0

def runSimple(n_levels, n_cells, output_filename='simple_run_output.txt'):

    # Write ForestClaw input file
    # contents = writeFileContents(
    #     clawpatch_mx=n_cells,
    #     clawpatch_my=n_cells,
    #     options_minlevel=n_levels,
    #     options_maxlevel=n_levels
    # )

    # os.remove('fclaw_options.ini')
    # with open('fclaw_options.ini', 'w') as file:
    #     file.write(contents)

    # Run executable with input file
    commands = [
        './simple',
        'fclaw_options.ini',
        '--report-timing-verbosity', 'all',
        '--clawpatch:mx', str(n_cells),
        '--clawpatch:my', str(n_cells),
        '--minlevel', str(n_levels),
        '--maxlevel', str(n_levels)
    ]
    with open(output_filename, 'w') as out_file:
        subprocess.run(commands, stdout=out_file)
    
    return

def extractResults(output_filename):

    results = SimpleResults()

    # Open output
    with open(output_filename, 'r') as in_file:

        # Read lines to get line number for useful info
        lines = in_file.readlines()
        for line_number, line in enumerate(lines):
            if '[fclaw] ----- End of HPS solver -----' in line:
                error_line_number = line_number + 1
            if '[libsc] Statistics for   EXTRA1' in line:
                setup_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA2' in line:
                build_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA3' in line:
                solve_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA4' in line:
                copy_time_line_number = line_number + 2
        
        # Extract results
        results.L1_error = float(lines[error_line_number].split()[3])
        results.L2_error = float(lines[error_line_number].split()[4])
        results.LI_error = float(lines[error_line_number].split()[5])
        results.setup_time = float(lines[setup_time_line_number].split()[5])
        results.build_time = float(lines[build_time_line_number].split()[5])
        results.solve_time = float(lines[solve_time_line_number].split()[5])
        results.copy_time = float(lines[copy_time_line_number].split()[5])

    return results

def writeResults(results, file):

    line = ''
    line += '{},'.format(results.n_levels)
    line += '{},'.format(results.n_cells)
    line += '{},'.format(results.dofs)
    line += '{},'.format(results.L1_error)
    line += '{},'.format(results.L2_error)
    line += '{},'.format(results.LI_error)
    line += '{},'.format(results.setup_time)
    line += '{},'.format(results.build_time)
    line += '{},'.format(results.solve_time)
    line += '{},'.format(results.copy_time)
    line += '\n'

    file.write(line)

def computeDOFs(n_levels, n_cells):
    return 2**(2*n_levels-2) * n_cells**2

def main():

    # Variables
    # number_of_levels = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    # number_of_cells = np.array([4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
    number_of_levels = np.arange(1, 3)
    number_of_cells = 2**np.arange(2, 6)
    output_filename = 'simple_run_output.txt'
    max_dofs = 2**22 # About 4 million DOFs

    N, M = np.meshgrid(number_of_levels, number_of_cells)
    DOFs = computeDOFs(N, M)
    total_runs = np.sum(DOFs <= max_dofs)
    current_run = 0

    column_names = [
        'n_levels',
        'n_cells',
        'dofs',
        'L1_error',
        'L2_error',
        'LI_error',
        'setup_time',
        'build_time',
        'solve_time',
        'copy_time'
    ]

    try:
        os.remove('timing_and_error.csv')
    except FileNotFoundError:
        pass
    for n in number_of_levels:
        for m in number_of_cells:
            dofs = computeDOFs(n, m)
            if dofs <= max_dofs:

                print('Running `simple` - Run {} of {}'.format(current_run, total_runs))
                print('  n = {}, m = {}, dofs = {}'.format(n, m, dofs))

                # Run simple exectutable
                runSimple(n, m, output_filename=output_filename)

                # Extract necessary data
                results = extractResults(output_filename)
                results.n_levels = n
                results.n_cells = m
                results.dofs = computeDOFs(n, m)
                
                # Write data to file
                with open('timing_and_error.csv', 'a') as file:
                    if current_run == 0:
                        [file.write('{},'.format(c)) for c in column_names]
                        file.write('\n')
                    writeResults(results, file)
                
                current_run += 1

if __name__ == '__main__':
    main()