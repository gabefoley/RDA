from werkzeug.datastructures import FileStorage
from wtforms import ValidationError, fields
from wtforms.validators import required
from wtforms.widgets import FileInput
from flask_wtf import FlaskForm
from wtforms import SubmitField, BooleanField, DateField, FileField, SelectMultipleField, SelectField, FieldList, \
    FormField, \
                                                                          TextField, \
    SelectMultipleField, \
    TextAreaField, StringField
from wtforms.widgets import ListWidget, CheckboxInput
from wtforms.validators import DataRequired, Email
from gettext import gettext


# Form for uploading files
class UploadForm(FlaskForm):
    upload_seqs = FileField(u'Upload your sequences (FASTA format)', validators=[DataRequired()])
    upload_submit = SubmitField("Upload sequences ")


class AnnotationForm(FlaskForm):
    '''
    Form to specify which annotation form should be displayed
    '''
    database = SelectField('Which database do you wish to retrieve annotations from?', choices=[('UniProt',
                                                                                                 'UniProt'), ('BRENDA',
                                                                                                 'BRENDA')])

class DownloadForm(FlaskForm):
    '''
    Form to trigger the downloading of whatever has been annotated
    '''
    download = SubmitField('Download annotation file')


class BrendaForm(FlaskForm):
    brenda_upload = FileField('Upload your BRENDA annotation file', validators=[DataRequired()])
    brenda_species = SelectField('Map based on organism name or UniProt ID', choices=[('Organism name',
                                                                                      'Organism name'),
                                                                                      ('UniProt ID', 'UniProt ID')])
    brenda_ubiquitous = BooleanField('Add annotations common to all')
    brenda_comments = BooleanField('Add annotations found in comments')

    brenda_select = SelectMultipleField('Which BRENDA annotations should be retrieved?',
                                        choices=[("APPLICATION", "APPLICATION"),
                                                 ("ACTIVATING_COMPOUND", "ACTIVATING_COMPOUND"), ("CLONED", "CLONED"),
                                                 ("COFACTOR", "COFACTOR"), ("ENGINEERING", "ENGINEERING"),
                                                 ("EXPRESSION", "EXPRESSION"),
                                                 ("GENERAL_STABILITY", "GENERAL_STABILITY"),
                                                 ("IC50_VALUE", "IC50_VALUE"),
                                                 ("INHIBITORS", "INHIBITORS"), ("KI_VALUE", "KI_VALUE"),
                                                 ("KM_VALUE", "KM_VALUE"),
                                                 ("LOCALIZATION", "LOCALIZATION"), ("METALS_IONS", "METALS_IONS"),
                                                 ("MOLECULAR_WEIGHT", "MOLECULAR_WEIGHT"),
                                                 ("NATURAL_SUBSTRATE_PRODUCT", "NATURAL_SUBSTRATE_PRODUCT"),
                                                 ("OXIDATION_STABILITY", "OXIDATION_STABILITY"),
                                                 ("PH_OPTIMUM", "PH_OPTIMUM"),
                                                 ("PH_RANGE", "PH_RANGE"), ("PURIFICATION", "PURIFICATION"),
                                                 ("REACTION", "REACTION"),
                                                 ("REACTION_TYPE", "REACTION_TYPE"),
                                                 ("RECOMMENDED_NAME", "RECOMMENDED_NAME"),
                                                 ("REFERENCE", "REFERENCE"), ("SOURCE_TISSUE", "SOURCE_TISSUE"),
                                                 ("SPECIFIC_ACTIVITY", "SPECIFIC_ACTIVITY"),
                                                 ("STORAGE_STABILITY", "STORAGE_STABILITY"),
                                                 ("SUBSTRATE_PRODUCT", "SUBSTRATE_PRODUCT"), ("SUBUNITS", "SUBUNITS"),
                                                 ("SYNONYMS", "SYNONYMS"), ("SYSTEMATIC_NAME", "SYSTEMATIC_NAME"),
                                                 ("TEMPERATURE_OPTIMUM", "TEMPERATURE_OPTIMUM"),
                                                 ("TEMPERATURE_RANGE", "TEMPERATURE_RANGE"),
                                                 ("TEMPERATURE_STABILITY", "TEMPERATURE_STABILITY"),
                                                 ("TURNOVER_NUMBER", "TURNOVER_NUMBER")], validators=[DataRequired()])

    brenda_retrieve = SubmitField("Retrieve BRENDA annotations")

class UniProtForm(FlaskForm):
    uniprot_select = SelectMultipleField('Which UniProt annotations should be retrieved?', choices=[
('entry name', 'Entry name'),
('genes', 'Gene names'),
('genes(PREFERRED)', 'Gene names (primary)'),
('genes(ALTERNATIVE)', 'Gene names (synonym)'),
('genes(OLN)', 'Gene names (ordered locus)'),
('genes(ORF)', 'Gene names (ORF)'),
('organism', 'Organism'),
('organism-id', 'Organism ID'),
('protein names', 'Protein names'),
('proteome', 'Proteomes'),
('lineage(ALL)', 'Taxonomic lineage'),
('virus hosts', 'Virus hosts'),
('fragment', 'Fragment'),
('encodedon', 'Gene encoded by'),
('comment(ALTERNATIVE PRODUCTS)', 'Alternative products'),
('comment(ERRONEOUS GENE MODEL PREDICTION)', 'Erroneous gene model prediction'),
('comment(ERRONEOUS INITIATION)', 'Erroneous initiation'),
('comment(ERRONEOUS TERMINATION)', 'Erroneous termination'),
('comment(ERRONEOUS TRANSLATION)', 'Erroneous translation'),
('comment(FRAMESHIFT)', 'Frameshift'),
('comment(MASS SPECTROMETRY)', 'Mass spectrometry'),
('comment(POLYMORPHISM)', 'Polymorphism'),
('comment(RNA EDITING)', 'RNA editing'),
('comment(SEQUENCE CAUTION)', 'Sequence caution'),
('length', 'Length'),
('mass', 'Mass'),
('sequence', 'Sequence'),
('feature(ALTERNATIVE SEQUENCE)', 'Alternative sequence'),
('feature(NATURAL VARIANT)', 'Natural variant'),
('feature(NON ADJACENT RESIDUES)', 'Non-adjacent residues'),
('feature(NON STANDARD RESIDUE)', 'Non-standard residue'),
('feature(NON TERMINAL RESIDUE)', 'Non-terminal residue'),
('feature(SEQUENCE CONFLICT)', 'Sequence conflict'),
('feature(SEQUENCE UNCERTAINTY)', 'Sequence uncertainty'),
('version(sequence)', 'Sequence version'),
('ec', 'EC number'),
('comment(ABSORPTION)', 'Absorption'),
('comment(CATALYTIC ACTIVITY)', 'Catalytic activity'),
('chebi', 'ChEBI'),
('chebi(Catalytic activity)', 'ChEBI (Catalytic activity)'),
('chebi(Cofactor)', 'ChEBI (Cofactor)'),
('chebi-id', 'ChEBI IDs'),
('comment(COFACTOR)', 'Cofactor'),
('comment(ENZYME REGULATION)', 'Enzyme regulation'),
('comment(FUNCTION)', 'Function[CC]'),
('comment(KINETICS)', 'Kinetics'),
('comment(PATHWAY)', 'Pathway'),
('comment(REDOX POTENTIAL)', 'Redox potential'),
('comment(TEMPERATURE DEPENDENCE)', 'Temperature dependence'),
('comment(PH DEPENDENCE)', 'pH dependence'),
('feature(ACTIVE SITE)', 'Active site'),
('feature(BINDING SITE)', 'Binding site'),
('feature(DNA BINDING)', 'DNA binding'),
('feature(METAL BINDING)', 'Metal binding'),
('feature(NP BIND)', 'Nucleotide binding'),
('feature(SITE)', 'Site'),
('annotation score', 'Annotation score'),
('features', 'Features'),
('comment(CAUTION)', 'Caution'),
('comment(MISCELLANEOUS)', 'Miscellaneous[CC]'),
('keywords', 'Keywords'),
('context', 'Matched text'),
('existence', 'Protein existence'),
('tools', 'Tools'),
('reviewed', 'Reviewed'),
('comment(SUBUNIT)', 'Subunit structure[CC]'),
('interactor', 'Interacts with'),
('comment(DEVELOPMENTAL STAGE)', 'Developmental stage'),
('comment(INDUCTION)', 'Induction'),
('comment(TISSUE SPECIFICITY)', 'Tissue specificity'),
('go', 'Gene ontology (GO)'),
('go(biological process)', 'Gene ontology (biological process)'),
('go(molecular function)', 'Gene ontology (molecular function)'),
('go(cellular component)', 'Gene ontology (cellular component)'),
('go-id', 'Gene ontology IDs'),
('comment(ALLERGEN)', 'Allergenic properties'),
('comment(BIOTECHNOLOGY)', 'Biotechnological use'),
('comment(DISRUPTION PHENOTYPE)', 'Disruption phenotype'),
('comment(DISEASE)', 'Involvement in disease'),
('comment(PHARMACEUTICAL)', 'Pharmaceutical use'),
('comment(TOXIC DOSE)', 'Toxic dose'),
('comment(SUBCELLULAR LOCATION)', 'Subcellular location[CC]'),
('feature(INTRAMEMBRANE)', 'Intramembrane'),
('feature(TOPOLOGICAL DOMAIN)', 'Topological domain'),
('feature(TRANSMEMBRANE)', 'Transmembrane'),
('comment(PTM)', 'Post-translational modification'),
('feature(CHAIN)', 'Chain'),
('feature(CROSS LINK)', 'Cross-link'),
('feature(DISULFIDE BOND)', 'Disulfide bond'),
('feature(GLYCOSYLATION)', 'Glycosylation'),
('feature(INITIATOR METHIONINE)', 'Initiator methionine'),
('feature(LIPIDATION)', 'Lipidation'),
('feature(MODIFIED RESIDUE)', 'Modified residue'),
('feature(PEPTIDE)', 'Peptide'),
('feature(PROPEPTIDE)', 'Propeptide'),
('feature(SIGNAL)', 'Signal peptide'),
('feature(TRANSIT)', 'Transit peptide'),
('3d', '3D'),
('feature(BETA STRAND)', 'Beta strand'),
('feature(HELIX)', 'Helix'),
('feature(TURN)', 'Turn'),
('citationmapping', 'Mapped PubMed ID'),
('citation', 'PubMed ID'),
('created', 'Date of creation'),
('last-modified', 'Date of last modification'),
('sequence-modified', 'Date of last sequence modification'),
('version(entry)', 'Entry version'),
('comment(DOMAIN)', 'Domain[CC]'),
('comment(SIMILARITY)', 'Sequence similarities'),
('families', 'Protein families'),
('feature(COILED COIL)', 'Coiled coil'),
('feature(COMPOSITIONAL BIAS)', 'Compositional bias'),
('feature(DOMAIN EXTENT)', 'Domain[FT]'),
('feature(MOTIF)', 'Motif'),
('feature(REGION)', 'Region'),
('feature(REPEAT)', 'Repeat'),
('feature(ZINC FINGER)', 'Zinc finger'),
('lineage(all)', 'Taxonomic lineage (all)'),
('lineage(SUPERKINGDOM)', 'Taxonomic lineage (SUPERKINGDOM)'),
('lineage(KINGDOM)', 'Taxonomic lineage (KINGDOM)'),
('lineage(SUBKINGDOM)', 'Taxonomic lineage (SUBKINGDOM)'),
('lineage(SUPERPHYLUM)', 'Taxonomic lineage (SUPERPHYLUM)'),
('lineage(PHYLUM)', 'Taxonomic lineage (PHYLUM)'),
('lineage(SUBPHYLUM)', 'Taxonomic lineage (SUBPHYLUM)'),
('lineage(SUPERCLASS)', 'Taxonomic lineage (SUPERCLASS)'),
('lineage(CLASS)', 'Taxonomic lineage (CLASS)'),
('lineage(SUBCLASS)', 'Taxonomic lineage (SUBCLASS)'),
('lineage(INFRACLASS)', 'Taxonomic lineage (INFRACLASS)'),
('lineage(SUPERORDER)', 'Taxonomic lineage (SUPERORDER)'),
('lineage(ORDER)', 'Taxonomic lineage (ORDER)'),
('lineage(SUBORDER)', 'Taxonomic lineage (SUBORDER)'),
('lineage(INFRAORDER)', 'Taxonomic lineage (INFRAORDER)'),
('lineage(PARVORDER)', 'Taxonomic lineage (PARVORDER)'),
('lineage(SUPERFAMILY)', 'Taxonomic lineage (SUPERFAMILY)'),
('lineage(FAMILY)', 'Taxonomic lineage (FAMILY)'),
('lineage(SUBFAMILY)', 'Taxonomic lineage (SUBFAMILY)'),
('lineage(TRIBE)', 'Taxonomic lineage (TRIBE)'),
('lineage(SUBTRIBE)', 'Taxonomic lineage (SUBTRIBE)'),
('lineage(GENUS)', 'Taxonomic lineage (GENUS)'),
('lineage(SUBGENUS)', 'Taxonomic lineage (SUBGENUS)'),
('lineage(SPECIES GROUP)', 'Taxonomic lineage (SPECIES GROUP)'),
('lineage(SPECIES SUBGROUP)', 'Taxonomic lineage (SPECIES SUBGROUP)'),
('lineage(SPECIES)', 'Taxonomic lineage (SPECIES)'),
('lineage(SUBSPECIES)', 'Taxonomic lineage (SUBSPECIES)'),
('lineage(VARIETAS)', 'Taxonomic lineage (VARIETAS)'),
('lineage(FORMA)', 'Taxonomic lineage (FORMA)')]
,
                                         validators=[DataRequired()])
    uniprot_retrieve = SubmitField("Retrieve UniProt annotations")


class GOForm(FlaskForm):
    go_term = FileField('Upload your BRENDA annotation file', validators=[DataRequired()])
    uniprot_retrieve = SubmitField("Retrieve UniProt annotations")
