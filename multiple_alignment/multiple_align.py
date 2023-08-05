import requests
import gget
import subprocess
import os
import urllib.request

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

# gene_names = ['CFTR','NCAM1']

def get_alphafold_download_link(uniprot_id):
    
    link_pattern = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'
        
    return link_pattern.format(uniprot_id)

def download_alphafold_prediction(uniprot_id):

    if len(uniprot_id) > 6 and len(uniprot_id) < 10: uniprot_id = uniprot_id[:6]
    
    url = get_alphafold_download_link(uniprot_id)
    # result = subprocess.run(['wget', url, '-O', 'result.pdb'])
    urllib.request.urlretrieve(url, 'result.pdb')
    # return result

def gene2ensg(gene_name, species='human'):
    base_url = f"https://rest.ensembl.org"
    endpoint = f"/xrefs/symbol/{species}/{gene_name}?"

    url = base_url + endpoint
    headers = {"Content-Type": "application/json"}

    response = requests.get(url, headers=headers)

    if response.ok:
        data = response.json()
        if data:
            ensembl_id = data[0]['id']
            return ensembl_id
        else:
            print(f"No Ensembl ID found for the gene '{gene_name}'.")
    else:
        print(f"Error: {response.status_code} - {response.reason}")

def get_uniprot_from_gget_out (gget_result):
    start = gget_result[0].find('uniprot_id') + 12
    end = gget_result[0].find('ensembl_id') - 1 
    return gget_result[0][start:end]

def visualize_protein(gene_names):

    uniprot_list = []
    genes = {}
    for gene_name in gene_names:
        print(gene_name)
        ensembl = gene2ensg(gene_name)
        seq = gget.seq(ensembl, translate=True)
        uniprot_list.append(get_uniprot_from_gget_out(seq))
        genes[gene_name] = seq[1]
    # Write protein sequences to a temporary FASTA file
    fasta_file = "protein_sequences.fasta"
    with open(fasta_file, "w") as f:
        for gene_name, sequence in genes.items():
            f.write(f">{gene_name}\n{sequence}\n")
    # Perform multiple sequence alignment using MUSCLE
    # output_file = "aligned_sequences.fasta"
    subprocess.run('muscle -align protein_sequences.fasta -output aligned_sequences.fasta', shell=True)
    os.system('python stable.py protein_sequences.fasta aligned_sequences.fasta > stable.fasta')

    # Read the aligned sequences from the output file
    aligned_sequences = SeqIO.parse('stable.fasta', "fasta")
    output_file2 = 'aligned_sequences2.txt'
    with open(output_file2, 'w') as out:
        for record in aligned_sequences:
            out.write('>' + str(record.id) + '\n')
            out.write(str(record.seq) + '\n')

    with open(output_file2, 'r') as file:
        data = file.readlines()

    seq_array = []
    for i in range(0,len(data),2):
        seq_array.append(data[i+1])

    total_sequences = len(uniprot_list)

    for itr, uni in enumerate(uniprot_list):
        
        download_alphafold_prediction(uni)

        dashes = 0
        data_new = []
        for i in range(len(data[2*itr+1])):
            if data[2*itr+1][i] == '-' or data[2*itr+1][i] == '\n': 
                dashes += 1
                continue
            
            mismatch_ct = 0
            for itr2, seq in enumerate(seq_array):
                if itr2 == itr: continue
                if seq[i] != data[2*itr+1][i]:
                    mismatch_ct += 1
                    
            # char = data[2][i]
            if mismatch_ct == total_sequences - 1:  
                colr = 0
                colg = 0
                colb = 0
            elif mismatch_ct > int(total_sequences/2):  
                colr = 255
                colg = 0
                colb = 0
            elif mismatch_ct == int(total_sequences/2): 
                colr = 255
                colg = 255
                colb = 0
            elif mismatch_ct < int(total_sequences/2): 
                colr = 0
                colg = 255
                colb = 0
            else:
                colr = 69
                colg = 69
                colb = 69
            data_new.append("{{start_residue_number: {}, end_residue_number: {}, color:{{r:{},g:{},b:{}}}}}".format(i-dashes,i-dashes,colr,colg,colb))

        data_new_collapse = '[' + ', '.join(data_new) + ']'

        with open ('final_html_template.html', 'r') as template_html:
            to_edit = template_html.readlines()

        to_edit[144] = data_new_collapse + '})}) \n'

        with open ('result.pdb', 'r') as pdb:
            pdb_array = pdb.readlines()

        for i in range(len(pdb_array)):
            to_edit.insert(205+i, pdb_array[i])

        with open('templates/' + gene_names[itr] + '.html', 'w') as output:
            for line in to_edit:
                output.write(line)        

