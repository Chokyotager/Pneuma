from data import Genome
import numpy as np

genome = Genome()

key = "CM004307.1"

nucleotide, annotation = genome.getContig(key, end_at=100000)

threshold = 1

clustering = list(genome.annotation_fields)
print(clustering)

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
    for contigs in contiguous_coordinates:
        contigs[2] = sum(contigs[2])/(contigs[1] - contigs[0] + 1)


    return contiguous_coordinates

annotation = np.swapaxes(annotation, 0, 1)
annotations = dict()

for i in range(len(clustering)):
    cluster = clustering[i]
    annotations[cluster] = computeNucleotideScores(annotation[i])

all_annotations = dict()
all_annotations[key] = annotations

gff_annotations = list()
id_numeral = 0
for sequence in all_annotations.keys():

    sequence_annotations = all_annotations[sequence]

    for cluster in clustering:

        annotations = sequence_annotations[cluster]

        for annotation in annotations:

            data = "\t".join([sequence, "Pneuma", cluster, str(annotation[0]), str(annotation[1]), str(annotation[2]), ".", ".", "ID={};".format(id_numeral)])
            gff_annotations.append(data)
            id_numeral += 1

open("test/ground_truth.gff", "w+").write("\n".join(gff_annotations))
