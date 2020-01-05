import tensorflow as tf
import math
import random

from config import config

from model import createModel
from data import Genome

genome = Genome()
model = createModel(output_nodes=len(genome.annotation_fields))

max_frames = config["genome"]["max-scaffolds-per-batch"]

# Run through genome
epoch = int()
while True:

    if epoch >= config["training"]["epoch-limit"]:
        break

    print("Running epoch #{}".format(str(epoch)))

    for contig in genome.contigs:

      print("################ Running framed contig: {} ################".format(contig))

      size = random.randint(config["genome"]["min-scaffold-length"], config["genome"]["max-scaffold-length"])
      frame_limit = math.floor(config["genome"]["frame-adjustment-ratio"] * size)

      frame = random.randint(0, frame_limit)

      # Set previous frame

      more_frames = True
      previous_frame = 0

      while more_frames:
        sequence, annotations, more_frames, previous_frame = genome.getContigFrame(contig=contig, size=size, frame=frame, start_at=previous_frame, max_frames=max_frames, return_frame=True)

        if more_frames:
          status = "truncated: {} bp (frame-exclusive)".format(previous_frame)

        else:
          status = "end of sequence"

        print("Fitting framed contig: {} ({}) [contigs: {}, size/contig: {}, frame: {}/{}]".format(contig, status, len(sequence), size, frame, frame_limit))
        history = model.fit(x=sequence, y=annotations)

      model.save_weights(config["files"]["save-to"])

    epoch += 1
