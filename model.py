import tensorflow as tf
from tensorflow.keras import layers

from config import config
training_params = config["training"]

def createModel (output_nodes=None, compile=True):

    model = tf.keras.Sequential()

    gru_size = [5, 2]
    ff_size = [32, 4]

    model.add(layers.Input(batch_shape=[None, None, 5]))

    for i in range(gru_size[1]):
      model.add(layers.GRU(gru_size[0], return_sequences=True, reset_after=True, activation="selu", recurrent_activation="tanh", stateful=False))

    for i in range(ff_size[1]):
      model.add(layers.Dense(ff_size[0], activation="selu"))

    model.add(layers.Dense(output_nodes, activation="sigmoid"))

    model.summary()

    optimiser = tf.keras.optimizers.Adam(learning_rate=training_params["learning-rate"], beta_1=training_params["beta_1"], beta_2=training_params["beta_2"], amsgrad=True)

    if compile:
        model.compile(loss="binary_crossentropy", optimizer=optimiser, metrics=["mse", "accuracy"])

    return model
