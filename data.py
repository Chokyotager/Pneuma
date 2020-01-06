import numpy as np
import random
from Bio import SeqIO

from config import config

available_inputs = ["A", "T", "C", "G", "N"]
nucleotide_one_hot = np.eye(len(available_inputs))

def nucleotideToOneHot (nucleotide):
    nucleotide = nucleotide.upper()

    if nucleotide not in available_inputs:
        nucleotide = "N"

    index = available_inputs.index(nucleotide)

    return list(nucleotide_one_hot[index])

class Genome ():

  def __init__ (self):

    # Parse genome into memory
    self.sequences = list(SeqIO.parse(config["files"]["sequences"], "fasta"))

    # Parse GFF3/GFF/GTF/GBFF (GFX) file into nested dict format
    annotations = [x.split("\t") for x in open(config["files"]["annotations"]).read().split("\n") if len(x.split("\t")) > 2]

    self.annotations = dict()

    contigs = set([x[0] for x in annotations])
    annotation_fields = set([x[2].lower() for x in annotations])

    self.contigs = contigs
    self.annotation_fields = sorted(annotation_fields)

    for contig in contigs:

      contig_dict = dict()
      for annotation_field in annotation_fields:
        contig_dict[annotation_field] = list()

      self.annotations[contig] = contig_dict

    for annotation in annotations:

      contig = annotation[0]
      annotation_field = annotation[2].lower()

      start_coordinate = int(annotation[3])
      end_coordinate = int(annotation[4])

      self.annotations[contig][annotation_field].append([start_coordinate, end_coordinate])

  def getContig (self, name=None, start_at=0, end_at=None):

    # If none, randomise
    if name == None:
      name = random.choice(list(self.contigs))

    found_target = False
    for sequence in self.sequences:
      if name == sequence.name:
        target_sequence = sequence
        found_target = True
        break

    if found_target == False:
      raise ValueError("Invalid target sequence/contig name")

    target_sequence = str(target_sequence.seq)
    annotation = self.annotations[name]

    # Iterate through numerically
    # Perform maximum time-batch learning

    nucleotide_vector = list()
    annotation_vector = list()

    if end_at != None:
      end_at = min(end_at, len(target_sequence))

    else:
      end_at = len(target_sequence)

    # Optimisation
    last_iterations = [0] * len(self.annotation_fields)

    # Iterate through nucleotide
    for i in range(start_at, end_at):

      nucleotide = target_sequence[i]

      one_hot = nucleotideToOneHot(nucleotide)
      nucleotide_vector.append(list(one_hot))

      coordinate = i + 1

      # Search annotations
      annotation = [0] * len(self.annotation_fields)
      annotation_fields = list(self.annotation_fields)

      annotations = self.annotations[name]

      for j in range(len(annotation_fields)):
        annotation_field = annotation_fields[j]
        for k in range(last_iterations[j], len(annotations[annotation_field])):
          positions = annotations[annotation_field][k]

          #print(positions)

          # If within
          if coordinate >= positions[0] and coordinate <= positions[1]:
            annotation[j] = 1
            last_iterations[j] = max(0, k - 1)
            break

          if coordinate > positions[1]:
            last_iterations[j] = max(0, k - 1)

      annotation_vector.append(annotation)

    return nucleotide_vector, annotation_vector

  def getContigFrame (self, contig=None, size=80, start_at=0, frame=0, max_frames=None, return_frame=False):

    # Get contig sequence and splice
    if max_frames != None:
      contig, annotations = self.getContig(contig, start_at=start_at, end_at=((max_frames + 1) * size) + start_at)

    else:
      contig, annotations = self.getContig(contig, start_at=start_at)

    contig_frames = list()
    annotation_frames = list()

    more_frames = False
    frame_value = 0

    for i in range(frame, len(contig), size):

      if max_frames == 0:
        more_frames = True
        frame_value = i - frame + start_at
        break

      if len(contig) < (i + size):

        # Pad with Ns
        delta = i + size - len(contig)

        contig = contig + [[0, 0, 0, 0, 1]] * delta
        annotations = annotations + [[0] * len(annotations[0])] * delta

      contig_frames.append(contig[i : i + size])
      annotation_frames.append(annotations[i : i + size])

      max_frames = max_frames - 1

    if return_frame:

      return contig_frames, annotation_frames, more_frames, frame_value

    else:

      return contig_frames, annotation_frames

class UnclassifiedGenome ():

    def __init__ (self, input_file):

        self.sequences = list(SeqIO.parse(input_file, "fasta"))
        self.contigs = [x.name for x in self.sequences]

    def getContig (self, contig=None, start_at=None, end_at=None):

        if contig == None:
            contig = random.choice(self.contigs)

        found_target = False

        for sequence in self.sequences:
            if sequence.name == contig:
                found_target = True
                break

        if found_target == False:
            raise ValueError("Invalid target sequence/contig name")

        sequence = str(sequence.seq)

        nucleotide_vector = list()

        if start_at == None:
            start_at = 0

        if end_at == None:
            end_at = len(sequence)

        else:
            end_at = min(len(sequence), end_at)

        for i in range(start_at, end_at):
            vector = nucleotideToOneHot(sequence[i])
            nucleotide_vector.append(vector)

        return nucleotide_vector
