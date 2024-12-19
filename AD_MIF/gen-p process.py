import os
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
import time

#AD Gene phenotype Processing

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


def train_autoencoder(data, input_dim, hidden_dim, epochs=100, batch_size=64, lr=1e-3, log_file=None):
    data = torch.FloatTensor(data)
    dataset = TensorDataset(data, data)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model = Autoencoder(input_dim=input_dim, hidden_dim=hidden_dim)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    if log_file:
        log_file.write(f"Training Log for Autoencoder on Data Chunk\n")
        log_file.write(f"Data Dimensions: Input dim={input_dim}, Hidden dim={hidden_dim}\n")
        log_file.write("=================\n")
        log_file.write(f"Training start time: {start_time}\n")

    for epoch in tqdm(range(epochs), desc="Training Epochs"):
        for batch_data, _ in dataloader:
            if torch.isnan(batch_data).any():
                batch_data[torch.isnan(batch_data)] = 0
            output = model(batch_data)
            loss = criterion(output, batch_data)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        if log_file:
            log_file.write(f'Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.6f}\n')
        print(f'Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.6f}\n')

    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    if log_file:
        log_file.write(f"Training end time: {end_time}\n")
    return model


def extract_features(model, data):
    data = torch.FloatTensor(data)
    with torch.no_grad():
        low_dim_data = model.encoder(data)
    return low_dim_data.cpu().numpy()


folder_path = '../'
low_dim_data_all = []

log_file_path = 'training_log.txt'
with open(log_file_path, 'w') as log_file:

    code_start_time = time.time()
    log_file.write("Code Execution Start Time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    for chrom_num in tqdm(range(1, 2), desc="Processing Chromosomes"):
        file_name = f'filtered_genotype.chr{chrom_num}.csv'
        file_path = os.path.join(folder_path, file_name)

        if os.path.exists(file_path):
            print(f"Data Name：{file_name}")

            log_file.write(f"\nProcessing {file_name}\n")
            chunk_iter = pd.read_csv(file_path, chunksize=CHUNK_SIZE)

            chunk_results = []
            for i, chunk in enumerate(chunk_iter):
                print(f"Processing part {i + 1} of chromosome {chrom_num}, total {i + 1}/11 parts")
                input_data = chunk.values

                input_dim = input_data.shape[1]
                hidden_dim = 200
                if i == 0:
                    print(f"The data dimension of {file_name} is: {input_dim}")
                    log_file.write(f"The data dimension of {file_name} is: {input_dim}\n")

                if i == 0:
                    autoencoder_model = train_autoencoder(input_data, input_dim, hidden_dim, epochs=100, log_file=log_file)

                low_dim_data = extract_features(autoencoder_model, input_data)
                chunk_results.append(low_dim_data)

                log_file.write(f"Shape of extracted low-dimensional features: {low_dim_data.shape}\n")

            chrom_low_dim_data = np.vstack(chunk_results)
            log_file.write(
                f"Low-dimensional data shape after dimensionality reduction for chromosome {chrom_num}: {chrom_low_dim_data.shape}\n")
            low_dim_data_all.append(chrom_low_dim_data)

            low_dim_file_path = f'low_dim_filtered_genotype.chr{chrom_num}.csv'
            pd.DataFrame(chrom_low_dim_data).to_csv(low_dim_file_path, index=False)
            log_file.write(f"Low-dimensional data for chromosome {chrom_num} has been saved to '{low_dim_file_path}'\n")
            print(f"Low-dimensional data for chromosome {chrom_num} has been saved to '{low_dim_file_path}'")

            code_media_time = time.time()
            total_media_time = code_media_time - code_start_time
            print(f"Chromosome {chrom_num} processing completed, shape: {chunk.shape}")
            log_file.write(f"Chromosome {chrom_num} processing completed, shape: {chunk.shape}\n")

            log_file.write("End time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

            log_file.write(
                f"\nTime taken for this chromosome: {total_media_time // 60:.0f} minutes, {total_media_time % 60:.2f} seconds\n")

            import gc
            del chunk_iter
            gc.collect()

    code_end_time = time.time()
    total_time = code_end_time - code_start_time
    print("All chromosomes have been processed")
    log_file.write(
        f"\nAll chromosomes processed, total time: {total_time // 60:.0f} minutes, {total_time % 60:.2f} seconds\n")

low_dim_final = np.concatenate(low_dim_data_all, axis=1)
print(f"Final data shape after concatenation: {low_dim_final.shape}")

final_df = pd.DataFrame(low_dim_final)
final_df.to_csv('low_dim_AD_features.csv', index=False)
print("Dimensionality-reduced feature data has been saved to 'low_dim_AD_features.csv'")
print(f"All information has been logged to: {log_file_path}")
