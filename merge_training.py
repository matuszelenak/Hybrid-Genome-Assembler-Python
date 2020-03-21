import numpy as np
from keras.models import Sequential
from keras.layers import Dense

dataset = np.loadtxt('read_data/ecoli.fasta_art_150b_15x.fq__read_data.json__merge_cases.csv', delimiter=',')
np.random.shuffle(dataset)

X = dataset[:, 0:5]
y = dataset[:, 5]

model = Sequential()
model.add(Dense(10, input_dim=5, activation='relu'))
model.add(Dense(20, activation='relu'))
model.add(Dense(1, activation='sigmoid'))

model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit(X, y, epochs=5, batch_size=10)
model.save('ecoli_model.hdf5')
_, accuracy = model.evaluate(X, y)
print('Accuracy: %.2f' % (accuracy * 100))

# print(model.predict(X).flatten().tolist())
