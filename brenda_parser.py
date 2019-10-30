from collections import defaultdict
import regex as re

class AnnotatedProtein():

    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.prot_annots = defaultdict(list)
        self.prot_annots["organism_split"] = " ".join(self.name.split(" ")[0:2])
        self.prot_annots["id"] = self.id

    def __str__(self):
        return self.name

    def add_annotation(self, field, annotation, prefix_text=""):

        if annotation.strip() != "":

            # annotation = "YES"

            self.prot_annots[field].append(prefix_text + annotation.strip())

def parse_brenda(filename, annot_dict = {}, uniprot_mapping=False, add_commonalities=True, add_comments=True):

    # Keep a dictionary mapping the BRENDA assigned number to the AnnotatedProtein class
    num_to_prot = {}

    chunk = ""
    match = False

    # Regex for extracting the protein name
    prot_name_regex = re.compile('PR\t#(?P<prot_num>\d+)#(?P<prot_name>[\w\s.-]*)(?=\<)')

    # Regex for extracting the protein ID (if it exists)
    prot_id_regex = re.compile('[\d\w._]*\s(?=UniProt|GenBank|SwissProt)')

    # Regex for matching the carriage return / new line
    end_pattern = "\r*\n"

    # Regex for extracting the primary annotation and associated numbers
    annot_regex = re.compile('(?<=\w{2}[\t])(#(?P<primary_nums>(\d+(,)?)+)#\s+)?(?P<primary_annot>[\w\W]*?(?=>\n))')

    # Regex for extracting the comment annotation and associated number
    comment_regex = re.compile('(?P<comment_nums>(\d+(,)?)+)#\s*(?P<comment_annot>[^\<]*)')

    for line in open(filename):

        # Check to see if this is a protein definition
        if  prot_name_regex.match(line):

            hit = prot_name_regex.match(line)

            # Extract the protein information
            prot_name = hit.group('prot_name').strip()
            prot_num = hit.group('prot_num').strip()
            prot_id = prot_id_regex.search(prot_name).group(0).strip() if prot_id_regex.search(prot_name) else ""

            # Create the AnnotatedProtein
            annotated_prot = AnnotatedProtein(prot_id, prot_name)

            # Map the BRENDA assigned number to the AnnotatedProtein
            num_to_prot[prot_num] = annotated_prot

        # Else check to see if this is an annotation
        else:
            if line.strip() in annot_dict.keys():
                annot_key = line.strip()
                match = True
                continue

            # End condition to stop and process the annotation
            elif re.match(end_pattern, line):
                match = False
                if len(chunk) > 4: # Skip blank lines and annotations with no information

                    annots = annot_regex.finditer(chunk)

                    # Add the primary annotation to all of the identified proteins
                    for annot in annots:

                        common_to_all = True

                        if (annot.groupdict()['primary_nums']):
                            common_to_all = False

                        primary_num_list = [x for x in num_to_prot.keys()] if common_to_all else [x for x in
                                                                                                  annot.groupdict()['primary_nums'].split(',')]

                        primary_annot = annot.groupdict()['primary_annot']

                        # Remove any \n \t characters
                        primary_annot = re.sub('\n', '', primary_annot)


                        if "(#" in primary_annot:
                            comments = True
                            comment_str = primary_annot.split("(#")[1]
                            primary_annot = primary_annot.split("(#")[0]

                        else:
                            comments = False
                            primary_annot = primary_annot.split("<")[0]



                        for primary_num in primary_num_list:
                            if (common_to_all and add_commonalities) or not common_to_all:
                                if primary_num in num_to_prot: # Remove after finalising!!!!!!!!
                                    num_to_prot[primary_num].add_annotation(annot_key, primary_annot)

                        # If there are comments, process them
                        if add_comments and comments:

                                comments = comment_regex.finditer(comment_str)
                                for comment in comments:
                                    # if comment.groupdict()['comment_nums'] in num_to_prot:

                                    comment_num_list = [x for x in comment.groupdict()['comment_nums'].split(',')]
                                    comment_annot = comment.groupdict()['comment_annot']

                                    comment_annot = re.sub('\n', '', comment_annot)
                                    comment_annot = re.sub('\t', ' ', comment_annot)


                                    for comment_num in comment_num_list:
                                        if comment_num in num_to_prot: # Remove after finalising!!!!!!
                                            num_to_prot[comment_num].add_annotation(annot_key, "COMMENT: " +
                                                                                    comment_annot)

                chunk = ""

                continue

            elif match:
                chunk += line

    return num_to_prot


# annot_dict = { "NAME" : "", "KM_VALUE": "KM", "EXPRESSION": "EXP", "SYSTEMATIC_NAME" : "SN", "SUBSTRATE_PRODUCT":
#     "SP", "SPECIFIC_ACTIVITY" : "SA", "SOURCE_TISSUE" : "ST", "KI_VALUE" : "KI", "REACTION": "RE", "REACTION_TYPE":
#     "RT"}
#
# num_to_prot = parse_brenda('files/final.txt', annot_dict, add_commonalities=False, add_comments=True)
#
# for enzyme in num_to_prot.values():
#     print (enzyme.name)
#     for k,v in enzyme.prot_annots.items():
#         print (k,v)
#     print()
