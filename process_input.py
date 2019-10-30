from werkzeug.utils import secure_filename
import os
from rda import app
import sequence
from pandas import DataFrame


UPLOAD_FOLDER = os.getcwd() + "/static/uploads/"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def save_files(uploads):
    """
    Take the FileStorage objects in the input dictionary, save them to the server, and then return the updated mapping
    that maps to their file location
    :param uploads:
    :return:
    """

    seq_dict = {}

    for path,file in uploads.items():

        filename = secure_filename(file.filename)

        if filename != "":
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)


            seq_dict['filepath'] = filepath

    return seq_dict

def save_brenda_annotations(uploads):

    for path, file in uploads.items():

        filename = secure_filename(file.filename)

        if filename != "":
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

    return filepath

def get_ids_from_fasta(filepath):
    fasta = sequence.readFastaFile(filepath)

    seq_names = [x.name for x in fasta]
    seq_info = [x.info for x in fasta]
    # seq_info = [" ".join(zip(x.name, x.info)) for x in fasta]
    print (seq_names)
    print (seq_info)
    seq_array = DataFrame(seq_names, seq_info, columns=['id', 'info'])
    seq_array.set_index('id')
    return seq_array

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]