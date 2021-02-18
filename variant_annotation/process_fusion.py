import pandas as pd
import argparse
import sys
import os
import re

def main():
        parser = argparse.ArgumentParser(description='Fusion processing')
        parser.add_argument('-i', '--input', default=None, type=str, required=True,
                help="Mapping file")
        parser.add_argument('-d', '--database', default=None, type=str,
                help="Variant annotation database")
        parser.add_argument('-f', '--infor', default=None, type=str,
                help="Variant annotation information file")
        parser.add_argument('-o', '--out_file_prefix', default=None, type=str,
                help="Output file prefix")

        args = parser.parse_args(sys.argv[1:len(sys.argv)])
        
        map_file = args.input
        database = args.database
        infor = args.infor
        out_file_pre = args.out_file_prefix
        
        fusion_summary_file = out_file_pre + "_original_fusion_summary.txt"
        fusion_summary = open(out_file_pre + "_original_fusion_summary.txt", "w")
        new_infor = open(out_file_pre + "_varInfor_with_fusion.txt", "w")
        new_database = open(out_file_pre + "_varInfor_with_fusion.fasta", "w")
        dat = pd.read_csv(map_file, sep="\t", header=0, low_memory=False)

        regex = re.compile('[^a-zA-Z]')
        count = 0
        fusion_summary.write("key\tvariant_ID\tgene1\tgene2\tbreakpoint1\tbreakpoint2\tsequence\tsample\tvalue\n")

        for i, row in dat.iterrows():
                sample = row["sample"]
                file_row = row["file_path"]
                key_list = []
                with open(file_row) as data:
                        line = data.readline().strip()
                        line = data.readline().strip()
                        while line:
                                gene1 = line.split("\t")[0]
                                gene2 = line.split("\t")[1]
                                breakpoint1 = line.split("\t")[4]
                                breakpoint2 = line.split("\t")[5]
                                confidence = line.split("\t")[14]
                                sequence = line.split("\t")[28]
                                if sequence != "." and sequence != "":
                                        count += 1
                                        sequence = regex.sub('', sequence)
                                        sequence = sequence.upper()
                                        key = breakpoint1 + "|" + breakpoint2 + "|" + gene1 + "|" + gene2 + "|" + sequence + "sample:" + sample
                                        if key not in key_list:
                                                key_list.append(key)
                                                fusion_summary.write(breakpoint1 + "|" + breakpoint2 + "|" + gene1 + "|" + gene2 + "|" + sequence + "\t" + str(count) + "\t" + gene1 +"\t" + gene2 + "\t" + breakpoint1 + "\t" + breakpoint2 + "\t" + sequence + "\t" + "sample:" + sample + "\t" + "1" +"\n")
                                line = data.readline().strip()
        
        fusion_summary.close()
        dat_fusion = pd.read_csv(fusion_summary_file, sep="\t", header=0, low_memory=False, dtype=str)
        dat_fusion = dat_fusion.pivot(index='key',columns='sample',values='value')
        dat_fusion = dat_fusion.fillna('0')
        dat_fusion.columns.name = None
        dat_fusion = dat_fusion.reset_index()
        
        with open(database) as data:
                line = data.readline().strip()
                while line:
                    new_database.write(line + "\n")
                    line = data.readline().strip()

        sample_list = []
        line_count = -1
        with open(infor) as data:
                line = data.readline().strip()
                line_split = line.split("\t")
                for one_item in line_split:
                        if "sample" in one_item:
                                sample_list.append(one_item)
                new_infor.write(line + "\n")
                line = data.readline().strip()
                while line:
                    new_infor.write(line + "\n")
                    line_count += 1
                    line = data.readline().strip()
        for i, row in dat_fusion.iterrows():
                line_count += 1
                keys = row['key'].split("|")
                gene1 = keys[2]
                gene2 = keys[3]
                breakpoint1 = keys[0]
                breakpoint2 = keys[1]
                sequence = keys[4]
                
                new_database.write(">VAR|" + str(line_count) + "\n")
                new_database.write(sequence + "\n")
                
                new_infor.write("VAR|" + str(line_count) + "\tfusion\t" + breakpoint1 + "\t" + breakpoint2 + "\tfusion\tfusion\tfusion\tfusion\t" + gene1 + "/" + gene2 + "\tfusion\tfusion\tfusion\tfusion\tfusion\tfusion\t" + sequence)
                for one_sample in sample_list:
                        new_infor.write("\t" + str(row[one_sample]))
                new_infor.write("\n")
        new_infor.close()
        new_database.close()

if __name__ == '__main__':
        main()
