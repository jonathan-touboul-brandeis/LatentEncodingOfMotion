import os
import numpy as np
import torch
import torchvision.transforms as T
import h5py
from tqdm import tqdm
		
# Define the dimensions for each category
cat_dims = {
    'IJXX': (330, 510),
    'IJA2': (330, 498),
    'IJB2': (348, 510)
}

# Define the transformation to convert images to tensors
transform = T.Compose([
    T.ToTensor(),
])

def get_data(rootdir, cat, SF_indices, x, y, window_size,augment=False):
    filepath = os.path.join(rootdir, str(cat))					
    print('Loading at x=%d,y=%d, size: %d Done!',x,y,window_size)
    TF = []
    OR = []
    SFs = []
    images = []
    labels_TF = torch.arange(9).unsqueeze(1)
    labels_OR = torch.arange(200).unsqueeze(1)

    for SF in tqdm(range(SF_indices)):
        SF = SF + 1
        arrays = {}
        file = 'dataSF' + str(SF) + '.mat'
        f = h5py.File(filepath + '/' + file)
        for k, v in f.items():
            arrays[k] = np.array(v)
        a = np.reshape(np.transpose(arrays['signalS'], (2, 1, 0)),
                       (9, 200, cat_dims[str(cat[:4])][1], cat_dims[str(cat[:4])][0]))
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                a_norm = a[i, j, :, :]
                cropped_img = a_norm[x:x+window_size, y:y+window_size]
                if np.mean(cropped_img == 0) < 0.3:  # Drop if more than 30% black pixels
                    images.append(transform(cropped_img))
                    TF.append(labels_TF[i])
                    OR.append(labels_OR[j])
                    SFs.append(SF)
                    if augment:
                        images.append(transform(cropped_img))
                        TF.append(labels_TF[i])
                        OR.append(labels_OR[j])
                        SFs.append(SF)

    data = dict()
    data['images'] = images
    data['OR'] = OR
    data['TF'] = TF
    data['SF'] = SFs
    save_path = os.path.join(rootdir, str(cat) + '_data.pt')
    torch.save(data, save_path)
    
    
    return 'Done'


def get_full_data(rootdir, cat, SF_indices):
    filepath = os.path.join(rootdir, str(cat))
    TF = []
    OR = []
    SFs = []
    images = []
    labels_TF = torch.arange(9).unsqueeze(1)
    labels_OR = torch.arange(200).unsqueeze(1)

    for SF in tqdm(range(SF_indices)):
        SF = SF + 1
        arrays = {}
        file = 'dataSF' + str(SF) + '.mat'
        f = h5py.File(filepath + '/' + file)
        for k, v in f.items():
            arrays[k] = np.array(v)
        a = np.reshape(np.transpose(arrays['signalS'], (2, 1, 0)),
                       (9, 200, cat_dims[str(cat[:4])][1], cat_dims[str(cat[:4])][0]))
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                a_norm = a[i, j, :, :]
                images.append(transform(a_norm))
                TF.append(labels_TF[i])
                OR.append(labels_OR[j])
                SFs.append(SF)
    
    data = dict()
    data['images'] = images
    data['OR'] = OR
    data['TF'] = TF
    data['SF'] = SFs
    save_path = os.path.join(rootdir, str(cat) + '_fulldata.pt')
    torch.save(data, save_path)
    
    
    return data


def DatasetInWindow(dataset, x, y, window_size):
    data = dict()
    data['OR'] = dataset['OR']
    data['TF'] = dataset['TF']
    data['SF'] = dataset['SF']
    images = []
    for i in range(len(dataset['images'])):
        images.append(dataset['images'][i][:,x:x+window_size,y:y+window_size])
    data['images'] = images
    return data

def DatasetInWindowRescaled(dataset, x, y, window_size):
    data = dict()
    data['OR'] = dataset['OR']
    data['TF'] = dataset['TF']
    data['SF'] = dataset['SF']
    images = []
    for i in range(len(dataset['images'])):
        im=dataset['images'][i][:,x:x+window_size,y:y+window_size]
        im=(im-im.min())/(im.max()-im.min())
        images.append(im)
    data['images'] = images
    return data

def DatasetInWindowRescaled(dataset, x, y, window_size):
    data = dict()
    data['OR'] = dataset['OR']
    data['TF'] = dataset['TF']
    data['SF'] = dataset['SF']
    images = []
    for i in range(len(dataset['images'])):
        im=dataset['images'][i][:,x:x+window_size,y:y+window_size]
        im=(im-im.min())/(im.max()-im.min())
        images.append(im)
    data['images'] = images
    return data

def rescale_data(dataset):
    data = dict()
    data['OR'] = dataset['OR']
    data['TF'] = dataset['TF']
    data['SF'] = dataset['SF']
    images = []
    for i in range(len(dataset['images'])):
        im=dataset['images'][i]
        im=(im-im.min())/(im.max()-im.min())
        images.append(im)
    data['images'] = images
    return data


def Difference(dataset):
    data =dataset.copy()
    differences = []
    
    for i in range(len(dataset['images'])):
        FoundRef=0
        j=0
        while FoundRef==0 and j<len(dataset['images']):
            if dataset['TF'][j].item()==0 and dataset['SF'][j]==dataset['SF'][i] and dataset['OR'][j].item()==dataset['OR'][i].item():
                Reference=dataset['images'][j]
                FoundRef=1
            else:
                j+=1
        diff=dataset['images'][i]-Reference
        differences.append(diff)
    data['images'] = differences
    return data


def DifferenceRef(dataset,Ref):
    data =dataset.copy()
    differences = []
    
    for i in range(len(dataset['images'])):
        FoundRef=0
        j=0
        while FoundRef==0 and j<len(dataset['images']):
            if dataset['TF'][j].item()==Ref and dataset['SF'][j]==dataset['SF'][i] and dataset['OR'][j].item()==dataset['OR'][i].item():
                Reference=dataset['images'][j]
                FoundRef=1
            else:
                j+=1
        diff=dataset['images'][i]-Reference
        differences.append(diff)
    data['images'] = differences
    return data
