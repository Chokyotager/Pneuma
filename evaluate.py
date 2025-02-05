import numpy as np
import tensorflow as tf
import math
import random

from config import config

from model import createModel
from data import UnclassifiedGenome

# Generate a GFF file from the model output
genome = UnclassifiedGenome(input_file=config["files"]["evaluate"])
model = createModel(output_nodes=14)

clustering = ['transcript', 'centromere', 'trna', 'exon', 'region', 'snorna', 'sequence_feature', 'ncrna', 'mrna', 'gene', 'snrna', 'rrna', 'cds', 'pseudogene']
clustering = sorted(clustering)

annotations = dict()

model.load_weights(config["files"]["save-to"])

threshold = config["evaluation"]["threshold"]
minimum_basepairs = config["evaluation"]["minimum-basepairs"]

def computeNucleotideScores (vector):

    # [start, stop, % score]
    contiguous_coordinates = list()
    latest_contig = None

    for i in range(len(vector)):

        if vector[i] >= threshold:
            # Count

            # Add new entry
            if latest_contig == None or i > (latest_contig + 1):
                # Put to start of array
                contiguous_coordinates = [[i + 1, i + 1, [vector[i]]]] + contiguous_coordinates
                latest_contig = i

            else:
                contiguous_coordinates[0][1] = i + 1
                contiguous_coordinates[0][2].append(vector[i])
                latest_contig = i

    # Normalise contiguous scores
    qc_coordinates = list()
    for contigs in contiguous_coordinates:
        contigs[2] = sum(contigs[2])/(contigs[1] - contigs[0] + 1)

        if (contigs[1] - contigs[0] + 1) >= minimum_basepairs:
            qc_coordinates.append(contigs)


    return qc_coordinates

all_annotations = dict()

# Evaluate
for contig in genome.contigs:

    print("Evaluating contig: {}".format(contig))

    sequence = genome.getContig(contig, end_at=100000)

    output_vector = model.predict([sequence])[0]
    output_vector = np.swapaxes(output_vector, 0, 1)

    annotations = dict()
    for i in range(len(clustering)):
        cluster = clustering[i]
        annotations[cluster] = computeNucleotideScores(output_vector[i])

    all_annotations[contig] = annotations

# Generate annotations

gff_annotations = list()
id_numeral = 0

all_data = list()

for sequence in all_annotations.keys():

    sequence_annotations = all_annotations[sequence]

    for cluster in clustering:

        annotations = sequence_annotations[cluster]

        for annotation in annotations:

            all_data.append([sequence, "Pneuma", cluster, str(annotation[0]), str(annotation[1]), str(annotation[2]), ".", ".", "ID={};".format(id_numeral)])
            id_numeral += 1

def sortNumerically (element):
    return [element[0], int(element[4])]

all_data.sort(key=sortNumerically)

for datapoint in all_data:

    data = "\t".join(datapoint)
    gff_annotations.append(data)

open("test/evaluated.gff", "w+").write("\n".join(gff_annotations))
