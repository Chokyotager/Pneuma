import numpy as np
from Bio import SeqIO

from config import config

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
    self.annotation_fields = annotation_fields

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

    # Change to one-hot vectors with independent binary predictions
    available_inputs = ["A", "T", "C", "G", "N"]
    nucleotide_one_hot = np.eye(len(available_inputs))

    target_sequence = str(target_sequence.seq)
    annotation = self.annotations[name]

    # Iterate through numerically
    # Perform maximum time-batch learning

    nucleotide_vector = list()
    annotation_vector = list()

    if end_at != None:
      end_at = min(end_at, len(target_sequence))

    for i in range(start_at, end_at):
      nucleotide = target_sequence[i].upper()

      if (nucleotide not in available_inputs):
        nucleotide = "N"

      index = available_inputs.index(nucleotide)

      nucleotide_vector.append(list(nucleotide_one_hot[index]))

      coordinate = i + 1

      # Search annotations
      annotation = [0] * len(self.annotation_fields)
      annotation_fields = list(self.annotation_fields)

      annotations = self.annotations

      for j in range(len(annotation_fields)):
        annotation_field = annotation_fields[j]
        for k in range(len(annotations[name][annotation_field])):

          positions = annotations[name][annotation_field][k]

          if coordinate >= positions[0] and coordinate <= positions[1]:
            annotation[j] = 1
            annotations[name][annotation_field] = annotations[name][annotation_field][k - 1:]
            break

          if coordinate > positions[1]:
            annotations[name][annotation_field] = annotations[name][annotation_field][k - 1:]
            break

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