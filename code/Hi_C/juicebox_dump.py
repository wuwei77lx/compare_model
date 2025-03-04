import argparse
import subprocess
from subprocess import check_call

def parseargs():
    parser = argparse.ArgumentParser(description='Download and dump HiC data')
    parser.add_argument('--hic_file', required=True, help="Path or url to .hic file.")
    parser.add_argument('--juicebox', required=True, default="", help="path to juicebox executable or java command invoking juicer_tools.jar. eg: 'java -jar juicer_tools.jar'")
    parser.add_argument('--resolution', default=5000, help="Resolution of HiC to download. In units of bp.")
    parser.add_argument('--outdir', default=".")
    parser.add_argument('--include_raw', action="store_true", help="Download raw matrix in addtion to VC")
    parser.add_argument('--chromosomes', default="all", help="comma delimited list of chromosomes to download")
    parser.add_argument('--skip_gzip', action="store_true", help="dont gzip hic files")

    return parser.parse_args()

def run_command(command, **args):
    print("Running command: " + command)
    return check_call(command, shell=True, **args)

def main(args):

    if args.chromosomes == "all":
        chromosomes = list(range(1,20)) + ['X']
    else:
        chromosomes = args.chromosomes.split(",")

    for chromosome in chromosomes:
        print("Starting chr" + str(chromosome) + " ... ")
        outdir = "{0}/chr{1}".format(args.outdir, chromosome)
        command = "mkdir -p " + outdir
        out = subprocess.getoutput(command)

        ## Download observed matrix with VC normalization
        command = args.juicebox + " dump observed VC {0} chr{1} chr{1} BP {3} {2}/chr{1}.VCobserved".format(args.hic_file, chromosome, outdir, args.resolution)
        print(command)
        out = subprocess.getoutput(command)
        if not args.skip_gzip: 
            run_command("gzip {0}/chr{1}.VCobserved".format(outdir, chromosome))

        ## Download observed matrix with SCALE normalization
        command = args.juicebox + " dump observed SCALE {0} chr{1} chr{1} BP {3} {2}/chr{1}.SCALEobserved".format(args.hic_file, chromosome, outdir, args.resolution)
        print(command)
        out = subprocess.getoutput(command)
        if not args.skip_gzip: 
            run_command("gzip {0}/chr{1}.SCALEobserved".format(outdir, chromosome)) 

        ## Download VC norm file
        command = args.juicebox + " dump norm VC {0} chr{1} BP {3} {2}/chr{1}.VCnorm".format(args.hic_file, chromosome, outdir, args.resolution)
        out = subprocess.getoutput(command)
        print(command)
        if not args.skip_gzip: 
            run_command("gzip {0}/chr{1}.VCnorm".format(outdir, chromosome))
        
        ## Download SCALE norm file
        command = args.juicebox + " dump norm SCALE {0} chr{1} BP {3} {2}/chr{1}.SCALEnorm".format(args.hic_file, chromosome, outdir, args.resolution)
        out = subprocess.getoutput(command)
        print(command)
        if not args.skip_gzip: 
            run_command("gzip {0}/chr{1}.SCALEnorm".format(outdir, chromosome))
        
        if args.include_raw:
            ## Download raw observed matrix
            command = args.juicebox + " dump observed NONE {0} chr{1} chr{1} BP {3} {2}/chr{1}.RAWobserved".format(args.hic_file, chromosome, outdir, args.resolution)
            print(command)
            out = subprocess.getoutput(command)
            if not args.skip_gzip:
                run_command("gzip {0}/chr{1}.RAWobserved".format(outdir, chromosome))

if __name__ == '__main__':
    args = parseargs()
    main(args)
