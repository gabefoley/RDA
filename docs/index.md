# Retrieval of Database Annotations

[Retreivel of Database Annotations](https://mkdocs.org) is part of [GRASP-suite](https://grasp.scmb.uq.edu.au).

It is designed to add annotations from online databases to sets of protein sequences, in order for these annotations to be used within the ancestral sequence reconstruction program [GRASP](https://grasp.scmb.uq.edu.au).

## Basic operations

### UniProt annotations
* Upload a FASTA file with your protein sequences of interest
* Select the [UniProt database fields](https://www.uniprot.org/help/uniprotkb_column_names) you wish to add to your annotation file
* Download your annotation file and upload it into your chose reconstruction in [GRASP](https://grasp.scmb.uq.edu.au).

### BRENDA annotations
* Go to [BRENDA database](https://www.brenda-enzymes.org) and search for the protein of interest.
* At the bottom of this page, scroll down to the very bottom and select 
    * Download
    * BRENDA Textfile
* Upload this BRENDA annotations file (this is what is parsed by RDA)

### Combining annotations
It is fine to add annotations from multiple databases, simply select retrieve each set of annotations before downloading the annotation file. For a more detailed explanation of how combining annotations works please see the [Detailed explanation of combining annotations](#detailed) section. 

## BRENDA options

|  Option | Description  | Default  |
|---|---|---|
|  **Upload your BRENDA annotation file** | Upload the file which you want RDA to parse  | _Required_  |
|  **Map based on species name, not UniProt ID** |  BRENDA annotations sometimes contain information at the organism level, but not tied to specific UniProt records. <br><br>Selecting this option annotates the associated organisms in your uploaded FASTA file. <br><br>These are based on the organisms the species map to in UniProt. |  _True_ |
|  **Add annotations common to all** | When an annotation is found that is common to all proteins in the BRENDA text file (_not necessarily all proteins in the uploaded FASTA file_) add it  | _True_   |
| **Add annotations found in comments**  | Add annotations found in comment fields. Note: These can often contain important information   | _True_   |

<a name="splitting"></a>
## Splitting the organism name 

If you are choosing to annotate from BRENDA you have the option of mapping from organism name, not UniProt ID. In these cases, RDA only selects the first two words of the organism name from the BRENDA annotation file _and_ from the UniProt organism field in a means to better associate them, as each database sometimes appends additional information to this field.

For example -

```
Zymomonas mobilis subsp. mobilis P0DJA2 UniProt 
```

Would become

```
Zymomonas mobilis
```

To assess how well this has worked for your particular case, RDA outputs the following columns -

|  organism_uniprot  | organism_brenda  |organism_split|
|---|---|---|

So that any that are mapping incorrectly can be manually removed.

Q97UB2	Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2) (Sulfolobus solfataricus)
V9S9Y9	Sulfolobus acidocaldarius SUSAZ
A1B3F3 	Paracoccus denitrificans

<a name="detailed"></a>
## Detailed explanation of combining annotations

Lets assume we have the following files

**FASTA file** - seqs.fasta that has the following headers -

```
>tr|V9S9Y9|V9S9Y9_9CREN 
>sp|Q97UB2|ILVD_SULSO 
>tr|A1B3F3|A1B3F3_PARDP
```
**BRENDA annotations file** that has the following proteins -
```
PROTEIN
PR    #1# Salmonella enterica subsp. enterica serovar Typhimurium   <10,12,16>
PR    #2# Paracoccus denitrificans   <6>
PR    #21# Corynebacterium glutamicum Q8NQZ9 UniProt <34,36>
PR    #22# Sulfolobus solfataricus Q97UB2 SwissProt <28,39>
```

And we want to retrieve the annotation for `SUBSTRATE_PRODUCT` from BRENDA and the annotation for `mass` from UniProt

#### 1. First let's try just retrieving annotations from BRENDA based on UniProt ID

RDA can automatically retrieve the UniProt IDs from the **BRENDA annotations file**, but we can see from that file the only UniProt / SwissProt IDs that are present in that file are `Q8NQZ9` and `Q97UB2`.

And we can see that from our **FASTA file** only `Q97UB2` is present there.

So, our final annotations file to download will look like -

|  id | SUBSTRATE_PRODUCT_brenda |
|---|---|
|	Q97UB2|	FIXME substrate product|

**Note:we automatically add a suffix indicating which database the field refers to, in order to differentiate between them**

#### 2. Now let's try retrieving annotations from BRENDA based on the organism name

RDA can also match up the organism names between a **BRENDA annotations file** and a **FASTA file**. For more information see the [Splitting the organism name](#splitting) section.

RDA first maps all of the UniProt IDs in the **FASTA file** to their organism names, then splits the first two words from this annotation, resulting in a file like this -

|  id | organism_uniprot | organism_split	|
|---|---|---|
|	Q97UB2|	Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2) (Sulfolobus solfataricus)| Saccharolobus solfataricus	|
|	V9S9Y9|	Sulfolobus acidocaldarius SUSAZ|Sulfolobus acidocaldarius	|
|A1B3F3	| Paracoccus denitrificans	|Paracoccus denitrificans	|

From the **BRENDA annotations file** we can see that the organisms (after splitting) present are `Salmonella enterica`, `Paracoccus denitrificans`, `Corynebacterium glutamicum`, and `Sulfolobus solfataricus` - meaning we will be able to map both `Paracoccus denitrificans` and `Sulfolobus solfataricus` to their BRENDA annotations.

So, our final annotations file to download will look like (scroll to see all columns) -

|  id | organism_uniprot  | organism_brenda	| organism_split| SUBSTRATE_PRODUCT_brenda |
|---|---|---|---|---|
|Q97UB2	|Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2) (Sulfolobus solfataricus)	|	Sulfolobus solfataricus Q97UB2 SwissProt| Sulfolobus solfataricus	| FIXME_substrate	|
|A1B3F3	|Paracoccus denitrificans	|Paracoccus denitrificans	|Paracoccus denitrificans	| FIXME_substrate	|

**Note: If there are multiple UniProt IDs that map to the same organism name, all of them will received the BRENDA annotations associated with that organism**

#### 3. Now let's try retrieving annotations from UniProt

Retrievals from UniProt always work just on UniProt ID. Remember that we're looking for the annotation for `mass` from UniProt (as opposed to the annotation for `SUBSTRATE_PRODUCT` from BRENDA).

Because it is based on UniProt ID, we will be able to map any sequences in the **FASTA file** with a UniProt / SwissProt ID.

So the final annotations file to download will look like - 

|  id | mass_uniprot |
|---|---|
|	V9S9Y9|	FIXME mass|
|	Q97UB2|	FIXME mass|
|	A1B3F3|	FIXME mass|

**Note: We don't add in organism name FIXME: SHOULD WE?**

#### 4. Now let's try combining retrievals from UniProt and BRENDA

To combine retrievals, we upload a single FASTA file, add annotations from UniProt and BRENDA (in any order) _and then_ download the annotation file.

If we choose to combine a UniProt retrieval with a BRENDA based on UniProt ID, it will end up looking like this 

|  id | mass_uniprot | SUBSTRATE_PRODUCT_brenda |
|---|---|---|
|	V9S9Y9|	FIXME mass|  |
|	Q97UB2|	FIXME mass| FIXME substrate
|	A1B3F3|	FIXME mass||

**Note: As long as there is at least one annotation for a protein, it will be included**

If we choose to combine a UniProt retrieval with a BRENDA based on organism name, it will end up looking like this 

|  id | organism_uniprot | organism_brenda | organism_split| mass_uniprot | SUBSTRATE_PRODUCT_brenda|
|---|---|---|---|---|---|
|	V9S9Y9|	Sulfolobus acidocaldarius SUSAZ | | |FIXME mass | |
|	Q97UB2|	Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2) (Sulfolobus solfataricus)	| Sulfolobus solfataricus Q97UB2 SwissProt| Saccharolobus solfataricus |FIXME mass | FIXME substrate |
|	A1B3F3|	Paracoccus denitrificans| Paracoccus denitrificans |Paracoccus denitrificans | FIXME mass | FIXME substrate |


## Editing .tsv files

The downloaded annotation file is a .tsv file, we can easily be edited to remove, merge, or add to columns in any programming language or tool such as Microsoft's Excel.

Manually creating or adding in additional columns and annotations will also work and allow them to be mapped and queried by [GRASP](https://grasp.scmb.uq.edu.au). 


