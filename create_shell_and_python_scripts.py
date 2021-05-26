import os


def create_script(script_name):
    scripts = f'/dorotheeh/zivlang/scripts/'
    with open(f'{scripts}bash/{script_name}.sh', 'w') as f:
        f.write(f'''#!/bin/bash
#PBS -q dorotheeh
#PBS -N {script_name}
#PBS -e /dorotheeh/zivlang/err_and_out_files/{script_name}.ER
#PBS -o /dorotheeh/zivlang/err_and_out_files/{script_name}.OU
#PBS -l nodes=compute-0-311:ppn=10,mem=40000000kb

#eval "$(conda shell.bash hook)"
#conda activate zivlang

python /dorotheeh/zivlang/scripts/{script_name}.py
''')
    if not os.path.exists(scripts + script_name + '.py'):
        with open(scripts + script_name + '.py', 'w') as write_script:
            write_script.write('#!/usr/bin/env python\n')
        

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('script_name', help='script name')
    args = parser.parse_args()
    create_script(args.script_name)