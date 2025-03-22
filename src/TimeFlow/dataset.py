import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader

np.random.seed(1)
torch.manual_seed(1)

class CytometryData(Dataset):

    """
    CytometryData class reads a CSV file with flow or mass cytometry data, converts
    the numeric data into PyTorch tensors and splits it into train, validation and test sets
    for density estimation. To compare Real NVP architectures, modify the percentages of the
    train/validation/test sets must be modified.
    """

    def __init__(self, data, mode='train'):

        # Convert input data to numpy array
        if isinstance(data, pd.DataFrame):
            data = data.values.astype(np.float32)
        elif isinstance(data, np.ndarray):
            data = data.astype(np.float32)
        else:
            raise ValueError("Data should be saved as a pandas DataFrame or a numpy array.")

        # Calculate the number of rows for each split
        total_rows = data.shape[0]
        train_end = int(0.7 * total_rows)
        val_end = train_end + int(0.3 * total_rows)

        # Split input data into train/val/test sets
        if mode == 'train':
            self.data = torch.tensor(data[:train_end, :], dtype=torch.float32)

        elif mode == 'val':
            self.data = torch.tensor(data[train_end:val_end, :], dtype=torch.float32)

        elif mode == 'test':
            self.data = torch.tensor(data[val_end:, :], dtype=torch.float32)

        else:
            raise ValueError(f"Invalid mode {mode}, choose from 'train', 'val', 'test'")

    def __len__(self):

        return len(self.data)

    def __getitem__(self, idx):

        return self.data[idx]
