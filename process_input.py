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
            fasta = sequence.readFastaFile(filepath)
            seq_names = [x.name for x in fasta]
            seq_array = DataFrame(seq_names, columns=['id'])
            seq_dict[path] = filename
            seq_dict['seq_array'] = seq_array.to_msgpack()

    return seq_dict

def save_brenda_annotations(uploads):

    for path, file in uploads.items():

        filename = secure_filename(file.filename)

        if filename != "":
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

    return filepath