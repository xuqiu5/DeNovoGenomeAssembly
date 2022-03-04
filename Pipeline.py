#usr/bin/python
"""
Created on Sun Feb 20 20:32:07 2022

@author: Saideep Narendrula, HaojunSong
"""
import os
import subprocess
import argparse

def QC(input_dir, isolates_list, fastp_dir):
    # input directory should be full path to input files
    # fastp directory should be full path to fastp dir
    # run fastp 
    for iso in isolates_list:
        input_path_1 = input_dir + "/" + iso + "_1.fq.gz"
        input_path_2 = input_dir + "/" + iso + "_2.fq.gz"
        output_path_1 = fastp_dir + "/" + iso  + "_1_fastP.fq.gz"
        output_path_2 = fastp_dir + "/" + iso  + "_2_fastP.fq.gz"

        os.system("fastp -i {} -I {} -o {} -O {} --correction --cut_front --cut_tail -M 30 --thread 8 -f 5 -q 20 -e 28".format(input_path_1, input_path_2, output_path_1, output_path_2))
    
    return fastp_dir + "/"

def spades(input_dir,isolates_list, output_path):
    # output path should be working directory and then a folder titled spades
    os.system("cd {}".format(output_path))
    # make the output path directory
    os.mkdir(output_path+'/careful')
    # make contig output path directory
    os.mkdir(output_path+'/contigs')
    # Run spades
    for iso in isolates_list:
        input_path_1 = input_dir + "/" + iso + "_1_fastP.fq.gz"
        input_path_2 = input_dir + "/" + iso + "_2_fastP.fq.gz"
        output_path_1 = output_path + "/careful/" + iso 
        cmd_string = "spades.py -k 21,33,55,77  --careful  --pe1-1 {} --pe1-2 {} -o {}".format(input_path_1, input_path_2, output_path_1)
        os.system(cmd_string)
        
    # Give contigs their name
    for iso in isolates_list:
        output_path_1 = output_path + "/careful/" + iso + "/contigs.fasta"
        output_path_2 = output_path + "/careful/" + iso + "/${i}_contigs.fasta"
        cmd_string = "mv {} {}".format(output_path_1,output_path_2)
        os.system(cmd_string)
        
    # Move them to one folder
    for iso in isolates_list:
        output_path_1 = output_path + "/careful/" + iso + "/${i}_contigs.fasta"
        output_path_2 = output_path + "/contig/" + iso         
        cmd_string = "cp {} {}".format(output_path_1,output_path_2)
    
    # return location of contigs
    return output_path + "/contig/"

def idba(input_dir,isolates_list, output_path):
    # output path should be working directory and then a folder titled idba
    os.chdir(output_path)
    # make the fasta path directory
    os.mkdir("{}/fasta".format(output_path))
    # make idba output path directory
    os.mkdir("{}/idba_output".format(output_path))
    # make contig output path directory
    os.mkdir ("{}/contigs".format(output_path))
    
    os.system("cd {}".format(output_path + "/fasta"))
    # Uncompress and create fasta files
    for iso in isolates_list:
        # take the gz files and uncompress them for fq2fa (which only takes uncompressed files....)
        input_path_1 = input_dir + "/" + iso + "_1_fastP.fq.gz"
        input_path_2 = input_dir + "/" + iso + "_2_fastP.fq.gz"
        
        output_path_1 = iso + "_1_fastP.fq"
        output_path_2 = iso + "_2_fastP.fq"
        
        os.system("gzip -dc {} > {}".format(input_path_1, output_path_1))
        os.system("gzip -dc {} > {}".format(input_path_2, output_path_2))
        
       
        # make the fastq into fasta files and interleave them (thats the input type idba-ud accepts)                                
        output_path_3 =  iso +"_12.fas"
        os.system("fq2fa --merge {} {} {}".format(output_path_1, output_path_2, output_path_3))
        
    os.system("cd ..")
    # Run IDBA-UD
    for iso in isolates_list:
        # run idba-ud
        input_path_1 = "{}/fasta/{}_12.fas".format(output_path, iso) 
        output_path_1 = "{}/idba_output/{}".format(output_path,iso) 
        os.system("idba_ud -r {} -o {} --mink 23 --maxk 143 --step 10 --min_contig 300".format(input_path_1, output_path_1))
        
        # move to seperate folder
        input_path_2 = output_path_1 + "/contig.fa"
        output_path_2 = "{}/contigs/{}_contig.fa".format(output_path,iso)
        os.system("cp {} {}".format(input_path_2, output_path_2))
        
    return output_path + "/contigs/"
        

def megahit(input_dir,isolates_list, output_path):
    os.system("cd {}".format(output_path))
    # make contig output path directory
    os.mkdir(output_path+'/contigs')
    # Separate forward and backward strands for 150 bp
    forward1 = ""
    backward1 = ""
    count1 = 0
    total_150 = 0
    files1 = os.listdir(input_dir)
    for sample in files1:
        if "_1_fastP.fq.gz" in sample and count1 == 0:
            forward1 = sample
            count1 += 1
        if "_2_fastP.fq.gz" in sample and count1 == 1:
            backward1 = sample
            count1 += 1
        if count1 == 2:
            output_name = "MEGA_" + forward1[0:7]
            #if os.path.exists("/home/groupc/files/genome_assembly/megahit/" + output_name) == True:
            #    os.system("rm -r /home/groupc/files/genome_assembly/megahit/" + output_name)
            os.system("megahit -1 " + input_dir + "/" + forward1 + " -2 " + input_dir +"/" + backward1 + " -o " + output_name)
            os.system("mv " + output_name + " " + output_path + "/contigs/")
            forward1 = ""
            backward1 = ""
            count1 = 0
            total_150 += 1
            
    return output_path + "/contigs/"
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='A comprehensive all in one pipeline for genome assembly',
        epilog='The QC have to perform first if you have the raw sequence data. And this pipeline is tuned for BIOL7210 dataset. if you have trimmed dataset and only want to perform, you need to change the QC dir to designated place')
    parser.add_argument(
        '-i', '--input', dest="filepath", required=True, help='input folder path that are fasta files')
    parser.add_argument(
        '-qi', '--QCinput', dest="QCfilepath", required=False, help='input folder path that are trimmed fasta files')
    parser.add_argument(
        '-o', '--output', dest="output", type = str, required=True, help='output folder name')
    parser.add_argument(
        '-q', '--QC', dest='QC', action = 'store_true', required=False, help='will perform QC analysis if select ')
    parser.add_argument(
        '-a', '--assebly', dest='assembly', nargs='+', required=True, help='list the assembly that plan to use[spades, idba, megahit]')
  
    args = parser.parse_args()

    #getting the present working directory to perform analysis and store results
    current_path = os.getcwd()
    #the input dir/and isolates list used as parameter for each function
    input_dir = args.filepath
    #get only the isolate lists without extension
    isolates_list = []
    for file in os.listdir(input_dir):
        isolates_list.append(os.path.splitext(file)[0])
    # if trimmed data
    QC_output = args.QCfilepath

    os.mkdir(current_path +'/' + args.output)
    output_folder = current_path +'/' + args.output
    if args.QC:
        print('************** QC running ************')
        os.mkdir(output_folder + '/QC')
        QC_output = output_folder + '/QC'
        QC(input_dir, isolates_list, QC_output )
        print('************** QC done ************')
    if 'spades' in args.assembly:
        print('*************SPADES running *************')
        os.mkdir(output_folder + '/spades')
        spades_output = output_folder + '/spades'
        spades(QC_output, isolates_list, spades_output)
        print('*************SPADES done ***************')
    if 'IDBA' in args.assembly:
        print("***********IDBA runnings ***************")
        os.mkdir(output_folder + '/IDBA')
        IDBA_output = output_folder + '/IDBA'
        idba(QC_output, isolates_list, IDBA_output)
        print('*************IDBA done ************')
    #if 'Mega' in args.assembly:
    #    print('**********MEGA_hit Running************')
    #    os.mkdir(output_folder + '/mega')
    #    mega_output = output_folder + '/mega'
    #    megahit(QC_output, isolates_list, mega_output)
    #    print("**********MEGA_hit Done************")
        
    print('All requested runs completed')
