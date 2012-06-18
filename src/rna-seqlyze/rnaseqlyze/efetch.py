from Bio import Entrez

import rnaseqlyze
Entrez.email = rnaseqlyze.entrez_email

nc_db = "nuccore"
gb_type = "gb"
gb_mode = "text"

def get_nc_id(accession):
    handle = Entrez.esearch(db=nc_db, term=accession + "[Accession]")
    record = Entrez.read(handle)
    id_list = record["IdList"]
    if len(id_list) != 1:
        raise Exception("unexpected reply from Entrez: id_list: %s" % id_list)
    return id_list[0]

def fetch_nc_gb(gb_id, out_file):
    handle = Entrez.efetch(db=nc_db, id=gb_id, rettype=gb_type, retmode=gb_mode)
    while True:
        data = handle.read(4096)
        if data == '':
            break
        out_file.write(data)
