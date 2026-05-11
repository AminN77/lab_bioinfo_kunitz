import os, urllib.request

with open("/workspace/data/kunitz_pdb_nr20.fasta") as f:
    ids = [l[1:].strip() for l in f if l.startswith(">")]

os.makedirs("single_chain", exist_ok=True)
for entry in ids:
    pdb_id, chain = entry.split("_")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    raw = urllib.request.urlopen(url).read().decode()
    out = []
    for line in raw.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if len(line) > 21 and line[21] == chain:
                out.append(line)
        elif line.startswith(("HEADER","TITLE","END")):
            out.append(line)
    with open(f"single_chain/{pdb_id}_{chain}.pdb","w") as g:
        g.write("\n".join(out) + "\nEND\n")
    print(f"{pdb_id}_{chain}: {sum(1 for l in out if l.startswith('ATOM'))} ATOM lines")

print("Done.")
