import argparse
import shutil
import sys
import csv
import re
import subprocess
from typing import List
from pathlib import Path

TMPDIR = '/tmp'


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download fastq files from ENA",
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=35), 
        epilog="""
        
        """
    )
    
    parser.add_argument(
        "sample_file",
        nargs="?",
        default=None,
        help="sample.tsv file with SRR ids"
    )

    return parser.parse_args()


def get_srr_ids(sample_file) -> list[str]:
    srr = re.compile('^[SED]RR[0-9]*$')

    ids = []
    
    with open(sample_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            id = row['sample']
            m = srr.match(id)
            if m:
                ids.append(id)
        
    return ids


def download_srr(project_dir, id) -> int:
    cmd = ['fastq-dl', '-o', TMPDIR, '-a', id]
    try:
        result = subprocess.run(cmd)
        if result.returncode == 0:
            metadata = Path(TMPDIR) / 'fastq-run-info.tsv'
            se_src = Path(TMPDIR) / (id + '.fastq.gz')
            p1_src = Path(TMPDIR) / (id + '_1.fastq.gz')
            p2_src = Path(TMPDIR) / (id + '_2.fastq.gz')
            
            if metadata.exists():
                shutil.move(metadata, project_dir + '/reads/' + id + '.medadata.tsv')

            if(p1_src.exists() and p2_src.exists()):
                shutil.move(p1_src, project_dir + '/reads/' + id + '_1.fq.gz')
                shutil.move(p2_src, project_dir + '/reads/' + id + '_2.fq.gz')

            if(se_src.exists()):
                shutil.move(se_src, project_dir + '/reads/' + id + '_s.fq.gz')
                
            return 0
        return 1
    except Exception as e:
        print("ERROR: " + str(e))
        return 1


def downloaded(project_dir, id) -> bool:
    se_file = Path(project_dir + '/reads/' + id +'_s.fq.gz')
    p1_file = Path(project_dir + '/reads/' + id +'_1.fq.gz')
    p2_file = Path(project_dir + '/reads/' + id +'_2.fq.gz')

    return se_file.exists() or (p1_file.exists() and p2_file.exists())


def main():
    """Main entry point."""
    args = parse_arguments()

    sample_file = args.sample_file

    project_dir = re.sub('/reference/.*', '', sample_file)

    srr_ids = get_srr_ids(sample_file)

    for id in srr_ids:
        if not downloaded(project_dir, id):
            download_srr(project_dir, id)

    
if __name__ == "__main__":
    main()
