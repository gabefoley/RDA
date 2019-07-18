from flask import Flask, render_template, request, session, redirect, url_for, send_file, session, make_response, flash
from flask_bootstrap import Bootstrap
from models import UploadForm, AnnotationForm, UniProtForm, BrendaForm, DownloadForm
import pandas as pd
import brenda_parser
import brenda_annotations
import webservice
from wtforms.validators import StopValidation

import process_input


app = Flask(__name__)
Bootstrap(app)
app.config.from_pyfile('configs/config.py')

@app.route("/", methods=['GET', 'POST'])
def index():

    upload_form = UploadForm()
    annotation_form = AnnotationForm()
    uniprot_form = UniProtForm()
    brenda_form = BrendaForm()
    download_form = DownloadForm()

    if request.method == "POST":
        # print (upload_form.validate())
        # print (annotation_form.validate())
        # print (brenda_form.validate())
        # print (uniprot_form.validate())

        if upload_form.upload_submit.data:

            saved_files = process_input.save_files(request.files)
            session["data_generated"] = saved_files
            session['seq_array'] = saved_files['seq_array']

            return render_template("annotation.html", form=annotation_form, uniprot_form=uniprot_form,
                                   brenda_form=brenda_form, download_form=download_form, current_file=saved_files[
                'upload_seqs'])


        elif uniprot_form.uniprot_retrieve.data:

            # Get list of columns we want to annotate
            columns = request.form.getlist("uniprot_select")

            # Retrieve them from UniProt
            seq_array = pd.read_msgpack(session['seq_array'])
            id_list = seq_array['id'].values.tolist()

            print ('seq array ')
            print (seq_array)
            seq_array.set_index('id')


            # Check if user is specifically asking for organism and if it should be retained
            if 'organism' in columns:
                session['keep_organism'] = True

            # Check if we need to explicitly map to species (i.e. we haven't already or aren't just about to)
            if not 'organism' in seq_array and 'organism' not in columns:
                columns.insert(0, 'organism')



            # Retrieve the UniProt annotations
            uniprot_dict = webservice.getUniProtDict(id_list, columns)
            uniprot_array = pd.DataFrame.from_dict(uniprot_dict, orient='index')

            # Split up the organism field (to try and maximise matches between this field and the BRENDA annotations)
            uniprot_array['organism_split'] = uniprot_array.apply(lambda row: " ".join(row.organism.split(" ")[0:2]), axis=1)

            print ('before rename')

            print (uniprot_array)
            uniprot_array.index.names = ['id']
            # uniprot_array.reindex(uniprot_array.index.rename(['ID']))

            # uniprot_array = uniprot_array.rename(columns={uniprot_array.columns[0]: "ID"})

            session['uniprot_array'] = uniprot_array.to_msgpack()

            print ('uniprot dict ')
            for k,v in uniprot_dict.items():
                print (k,v)
            print()
            print ('uniprot array')
            print (uniprot_array)

            # Add the species information to the seq array
            if not 'organism' in seq_array:
                seq_array = pd.DataFrame.from_dict(uniprot_dict, orient='index', columns = ['organism',
                                                                                            'organism_split'] )

                print ('problem area is here ')
                print (seq_array)

                # seq_array = seq_array.rename(columns={seq_array.columns[0]: "id"})
                seq_array.index.name = 'id'

                # seq_array.set_index('id')

                print ('problem was renamed ')
                print (seq_array)


            session['seq_array'] = seq_array.to_msgpack()

            flash(u'Retrieved UniProt annotations', 'success')

            return render_template("annotation.html", form=annotation_form, uniprot_form=uniprot_form,
                                   brenda_form=brenda_form, download_form=download_form,
                                   selected='uniprot', current_file=session['data_generated']['upload_seqs'])

        elif brenda_form.brenda_retrieve.data:

            # List of sequences to add annotations to
            seq_array = pd.read_msgpack(session['seq_array'])


            # seq_array.index.names = ['id']

            print('here is seq array ')
            print(seq_array)

            # seq_array.set_index('id')


            # id_list = seq_array['id'].values.tolist()
            # id_list = seq_array['id'].values.tolist()

            id_list = seq_array.index.array


            # Get list of columns we want to get annotations for
            columns = request.form.getlist("brenda_select")

            # Get dictionary mapping annotations to abbreviated form
            annot_dict = brenda_annotations.get_brenda_dict(columns)

            # Get user's preferences for parsing the file
            brenda_species = True if 'brenda_species' in request.form else False
            session['brenda_species'] = brenda_species
            add_commonalities = True if 'add_commonalities' in request.form else False
            add_comments = True if 'add_comments' in request.form else False


            # Load BRENDA annotations file
            brenda_path = process_input.save_brenda_annotations(request.files)
            num_to_prot = brenda_parser.parse_brenda(brenda_path, annot_dict, brenda_species, add_commonalities,
                                                       add_comments)

            brenda_dict = {}

            for prot in num_to_prot.values():
                print (prot.name)
                print (prot.prot_annots)

                # If we're mapping based on species, add everything
                if brenda_species:
                    brenda_dict[prot.name] = prot.prot_annots

                # Else we need to check that our annotated ID from the BRENDA file is in our FASTA file
                else:
                    if prot.id in id_list:
                        brenda_dict[prot.name] = prot.prot_annots

            print ('brenda dict is ')
            for k,v in brenda_dict.items():
                print (k,v)

            brenda_array = pd.DataFrame.from_dict(brenda_dict, orient='index')

            # If mapping based on species name
            if brenda_species:

                # Check if we need to annotate the IDs from the FASTA with species names
                if not 'organism' in seq_array:
                    print ('Organism was not in seq array')
                    uniprot_dict = webservice.getUniProtDict(id_list, ['organism'])
                    for k, v in uniprot_dict.items():
                        print (k, v)
                    seq_array = pd.DataFrame.from_dict(uniprot_dict, orient='index', columns=['organism'])
                    print ('seq array before rename')
                    print (seq_array)
                    seq_array.index.names = ['id']

                    # seq_array = seq_array.rename(columns={seq_array.columns[0]: "ID"})
                    # seq_array['organism'] = 'str' + " ".join(seq_array['organism'].astype(str).split(" ")[0:2])

                    print ('Seq array is now ')
                    print (seq_array)
                    session['seq_array'] = seq_array.to_msgpack()

                else:
                    print ('Organism WAS IN seq array')

                # Merge the species and the BRENDA annotations
                brenda_array = pd.merge(seq_array, brenda_array, on=['organism'])

                # brenda_array = pd.merge(seq_array, brenda_array, on='organism', how='outer')

                print ('merged array is now ')
                print (brenda_array)


            # Otherwise we're mapping based on ID
            else:
                print ('Mapping based on ID')

                print ('seq array')
                print (seq_array)

                print (brenda_array)
                print (brenda_array)
                # Merge the FASTA IDs and the BRENDA annotations
                brenda_array = pd.merge(seq_array, brenda_array, on='id', how='outer')

                print ('BRENDA array is ')
                print (brenda_array)

                # Drop any entries without an ID
                brenda_array = brenda_array[pd.notnull(brenda_array['id'])]

                print ('And after dropping entries without an ID ')
                print (brenda_array)

            brenda_array.set_index('id')

            print ('brenda array before packing is ')
            print (brenda_array)

            session['brenda_array'] = brenda_array.to_msgpack()

            flash('Retrieved BRENDA annotations', 'success')

            return render_template("annotation.html", form=annotation_form, uniprot_form=uniprot_form,
                                   brenda_form=brenda_form, download_form=download_form,
                                   selected='brenda',
                                   current_file=session['data_generated'][
                'upload_seqs'])

        # Downloading the annotations file
        elif download_form.download.data:
            print ('Download the data')

            if 'uniprot_array' in session:
                print ('uniprot array was in session')
                uniprot_array = pd.read_msgpack(session['uniprot_array'])
                uniprot_array.columns = uniprot_array.columns.map(lambda x: x + '_uniprot' if x != 'id' and x != 'Name'
                else x)

            else:
                uniprot_array = pd.DataFrame()

            if 'brenda_array' in session:
                print ('brenda array was in session')
                print ('brenda after unpacking is ')
                brenda_array = pd.read_msgpack(session['brenda_array'])
                print (brenda_array)

                brenda_array.columns = brenda_array.columns.map(lambda x: x + '_brenda' if x != 'id' else x)
            else:
                brenda_array = pd.DataFrame()

            if uniprot_array.empty:
                if brenda_array.empty:
                    flash(u'Annotation file ended up being empty', 'error')
                    return render_template("upload.html", form=upload_form)
                else:
                    final_array = brenda_array
            else:
                if brenda_array.empty:
                    final_array = uniprot_array
                else:
                    if session['brenda_species']:
                        final_array = pd.merge(brenda_array, uniprot_array, on='organism_split', how='outer')
                    else:
                        final_array = pd.merge(brenda_array, uniprot_array, on='id', how='outer')


            print ('Final array is ')
            print (final_array)

            droplist = []
            dropcheck = ['organism_split_uniprot','organism_brenda', 'organism_split_brenda', 'organisim_x_brenda',
                         'organism_y_brenda', 'organism_x_uniprot', 'organism_y_uniprot']


            if 'brenda_species' in session:
                print ('brenda species in session')
                if not session['brenda_species']:
                    droplist = [i for i in final_array.columns if i in dropcheck]
            else:
                droplist = [i for i in final_array.columns if i in dropcheck]

            if 'keep_organism' in session:
                print ('keep organism in session')
                if not session['keep_organism']:
                    droplist += (i for i in ['organism_uniprot', 'organism_split_uniprot'] if i not in droplist and i
                                 in final_array.columns)
            else:
                droplist += (i for i in ['organism_uniprot', 'organism_split_uniprot'] if i not in droplist and i in
                             final_array.columns)

            print (droplist)
            print (type(droplist))

            final_array.drop(droplist, axis=1, inplace=True)

            # Drop any records that don't have at least 2 columns annotated (ID will always be annotated)
            # final_array = final_array.dropna(thresh=2)

            resp = make_response(final_array.to_csv(sep='\t', index=True))
            resp.headers["Content-Disposition"] = "attachment; filename=export.tsv"
            resp.headers["Content-Type"] = "text/tsv"

            return resp


    elif request.method == "GET":

        # Clear out any old annotations from previous sessions (good idea!)
        session.clear()

        return render_template("upload.html", form=upload_form)



@app.route("/about")
def about():
    return render_template("about.html")

if __name__ == "__main__":
    app.run(debug=True)