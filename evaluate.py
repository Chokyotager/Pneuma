import tensorflow as tf
import math
import random

from config import config

from model import createModel
from data import UnclassifiedGenome

# Generate a GFF file from the model output
genome = UnclassifiedGenome(input_file="data/GCF_000328475.2_Umaydis521_2.0_genomic.fna")
model = createModel(output_nodes=9)

clustering = list()

model.load_weights(config["files"]["save-to"])

# Evaluate
for contig in genome.contigs:

    contig = "NW_011929455.1"

    sequence = genome.getContig(contig)

    print(sequence)

    output_vector = model.predict([sequence])

    print(output_vector)

    break

# Evaluate accordingly
# GFF: [Sequence, source, feature, start, end, score, strand, phase, attributes]
