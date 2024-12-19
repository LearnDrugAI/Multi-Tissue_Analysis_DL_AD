import os
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm

#AD Gene exp Processing

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

CHUNK_SIZE = 50

class Autoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 1000),
            nn.ReLU(True),
            nn.Linear(1000, 500),
            nn.ReLU(True),
            nn.Linear(500, hidden_dim)
        )

        self.decoder = nn.Sequential(
            nn.Linear(hidden_dim, 500),
            nn.ReLU(True),
            nn.Linear(500, 1000),
            nn.ReLU(True),
            nn.Linear(1000, input_dim),
            nn.Sigmoid()
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded


def train_autoencoder(data, input_dim, hidden_dim, epochs=100, batch_size=64, lr=1e-3):
    data = torch.FloatTensor(data).to(device)
    dataset = TensorDataset(data, data)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model = Autoencoder(input_dim=input_dim, hidden_dim=hidden_dim).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    for epoch in tqdm(range(epochs), desc="Training Epochs"):
        for batch_data, _ in dataloader:
            if torch.isnan(batch_data).any():
                batch_data[torch.isnan(batch_data)] = 0
            output = model(batch_data)
            loss = criterion(output, batch_data)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        print(f'Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.6f}\n')

    return model


def extract_features(model, data):
    data = torch.FloatTensor(data).to(device)
    with torch.no_grad():
        low_dim_data = model.encoder(data)
    return low_dim_data.cpu().numpy()


gen_expression = np.load('gen_exp.npy')

input_dim = gen_expression.shape[1]
hidden_dim = 4400
print(f"gen_expression dim: {input_dim}")

autoencoder_model = train_autoencoder(gen_expression, input_dim, hidden_dim, epochs=2000)
low_dim_data = extract_features(autoencoder_model, gen_expression)

np.save('gen_exp_low_dim.npy', low_dim_data)
print("data saved as 'gen_exp_low_dim.npy'")
