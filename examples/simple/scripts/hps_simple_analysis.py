
# Imports
import sys
import shutil
import os
import subprocess

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

class SimpleResults:
    def __init__(self):
        self.solver = 'HPS'
        self.L1Error = 0
        self.L2Error = 0
        self.LIError = 0
        self.setupTime = 0
        self.buildTime = 0
        self.upwardsTime = 0
        self.solveTime = 0
        self.ellipticTime = 0
        self.nLevels = 0
        self.nCells = 0

def runSimple(nCellsSide, minLevel, maxLevel, iExample, outputFilename, runFISHPACK='F'):
    commands = [
        './simple',
        'fclaw_options.ini',
        '--report-timing', 'T',
        '--report-timing-verbosity', 'all',
        '--minlevel', str(minLevel),
        '--maxlevel', str(maxLevel),
        '--clawpatch:mx', str(nCellsSide),
        '--clawpatch:my', str(nCellsSide),
        '--clawpatch:refinement-criteria=value',
        '--user:example', str(iExample),
        '--hps:only_patch_solver={}'.format(runFISHPACK)
    ]
    with open(outputFilename, 'w') as out_file:
        subprocess.run(commands, stdout=out_file)

    return

def extractResults(outputFilename):

    results = SimpleResults()

    # Open output
    with open(outputFilename, 'r') as in_file:

        # Read lines to get line number for useful info
        lines = in_file.readlines()
        for line_number, line in enumerate(lines):
            if '[fclaw] error[0] =' in line:
                error_line_number = line_number
            if '[libsc] Statistics for   ELLIPTIC_SOLVE' in line:
                elliptic_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA1' in line:
                setup_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA2' in line:
                build_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA3' in line:
                upwards_time_line_number = line_number + 2
            if '[libsc] Statistics for   EXTRA4' in line:
                solve_time_line_number = line_number + 2
        
        # Extract results
        results.L1Error = float(lines[error_line_number].split()[3])
        results.L2Error = float(lines[error_line_number].split()[4])
        results.LIError = float(lines[error_line_number].split()[5])
        results.setupTime = float(lines[setup_time_line_number].split()[5])
        results.buildTime = float(lines[build_time_line_number].split()[5])
        results.upwardsTime = float(lines[upwards_time_line_number].split()[5])
        results.solveTime = float(lines[solve_time_line_number].split()[5])
        results.ellipticTime = float(lines[elliptic_line_number].split()[5])

    return results

def writeResults(results, file):

    line = ''
    line += '{},'.format(results.solver)
    line += '{},'.format(results.nLevels)
    line += '{},'.format(results.nCells)
    line += '{},'.format(results.L1Error)
    line += '{},'.format(results.L2Error)
    line += '{},'.format(results.LIError)
    line += '{},'.format(results.setupTime)
    line += '{},'.format(results.buildTime)
    line += '{},'.format(results.upwardsTime)
    line += '{},'.format(results.solveTime)
    line += '{},'.format(results.ellipticTime)
    line += '\n'

    file.write(line)

# Main function
def main(args):
    
    nCellsSideArray = 2**np.arange(4, 9)
    nLevelsRefinementArray = np.arange(0, 6)
    iExample = 4

    nRuns = len(nCellsSideArray) * len(nLevelsRefinementArray)

    resultsFileName = 'adaptive_hps_timing_and_error.csv'
    columnNames = [
        'solver',
        'nLevels',
        'nCells',
        'L1Error',
        'L2Error',
        'LIError',
        'setupTime',
        'buildTime',
        'upwardsTime',
        'solveTime',
        'ellipticTime'
    ]

    try:
        os.remove(resultsFileName)
    except FileNotFoundError:
        pass

    print('Begin HPS runs...')
    runCounter = 1
    for n in nCellsSideArray:
        for l in nLevelsRefinementArray:

            print('Running `simple` - run {} of {}: n = {}, l = {}'.format(runCounter, nRuns, n, l))

            # Run simple
            outputFileName = 'simple_n{}_l{}_e{}.txt'.format(n, l, iExample)
            runSimple(n, 1, l, iExample, outputFileName)

            # Extract necessary data
            results = extractResults(outputFileName)
            results.nCells = n
            results.nLevels = l
            results.solver = 'HPS'

            with open(resultsFileName, 'a') as file:
                if runCounter == 1:
                    [file.write('{},'.format(c)) for c in columnNames]
                    file.write('\n')
                writeResults(results, file)
            
            runCounter += 1
    print('Done!')

    print('Begin FISHPACK runs...')
    runCounter = 1
    for n in nCellsSideArray:
        for l in nLevelsRefinementArray:

            print('Running `FISHPACK` - run {} of {}: n = {}, l = {}'.format(runCounter, nRuns, n, l))

            # Run simple
            nCellsFISHPACK = n * 2**(l)
            outputFileName = 'FISHPACK_n{}_l{}_e{}.txt'.format(n, l, iExample)
            runSimple(nCellsFISHPACK, 0, 0, iExample, outputFileName, runFISHPACK='T')

            # Extract necessary data
            results = extractResults(outputFileName)
            results.nCells = nCellsFISHPACK
            results.nLevels = 1
            results.solver = 'FISHPACK'

            with open(resultsFileName, 'a') as file:
                writeResults(results, file)
            
            runCounter += 1
    print('Done!')

# Run main
if __name__ == '__main__':
    
    main(sys.argv)