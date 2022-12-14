# brain-topology

Code for the paper [Tree Representations of Brain Structural Connectivity via Persistent Homology](http://arxiv.org/abs/2211.00706.).

**Highlights of paper**:

* A novel application of computational topology to neuroimaging data.
* A new way to apply persistent homology ideas to brain connectome data.
* A new tree representation of brain networks.
* This representation reduces dimension relative to usual matrix representations.
* There are also clear interpretability advantages.
* Applications to HCP data leading to new insights into connectome-trait relationships.
* Opens up a new way to study brain structural connectomes in relation to other factors.
* Code are provided for easy implementation.

## Details

Below are the descriptions of folders in this repository:

* Folder `tools` includes code to convert a matrix of brain connections to our tree representation based on the Desikan-Killiany protocol and plotting fuctions to create a chordplot unique to our representation.

* Folder `data` includes data files used for tree construction and the Human Connectome Project's diffusion and structural MRI data. Some of the traits data in our analyses are restrictive. Interested readers may refer to [ConnectomeDB](https://wiki.humanconnectome.org/display/PublicData/How+to+Access+Data+on+ConnectomeDB) to gain access to download these files.

* Folder `paper` includes code to reproduce figures and analyses in our paper.