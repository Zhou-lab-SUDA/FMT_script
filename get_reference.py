import ete3, json
import sys
import click

def is_reference_genome(name):
    """Check if the leaf name is a reference genome (GCF_ followed by digits)"""
    if name.startswith("GCF_"):
        suffix = name[4:].replace('.', '')  # Remove any version dots
        return suffix.isdigit()
    return False


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
def main(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            if not parts or len(parts) < 5:
                continue
            
            tree_path = parts[0]
            sig = parts[1]
            rates = json.loads(parts[5])
            node_name = parts[4]
            
            try:
                # Load tree with flexible format detection
                tree = ete3.Tree(tree_path, format=1)
            except Exception as e:
                print(f"Error loading tree {tree_path}: {str(e)}", file=sys.stderr)
                fout.write(f"{tree_path}\t{node_name}\tNA\tNA\n")
                continue
                
            # Find the target node
            node = None
            for n in tree.traverse():
                if n.name == node_name:
                    node = n
                    break
            
            if node is None:
                print(f"Node {node_name} not found in {tree_path}", file=sys.stderr)
                fout.write(f"{tree_path}\t{node_name}\tNA\tNA\n")
                continue
            
            ingroup = node.get_leaf_names()
            in_ref = [g for g in ingroup if g.find('|') < 0]
            outgroup = set(tree.get_leaf_names()) - set(ingroup)
            out_ref = [g for g in outgroup if g.find('|') < 0]
            
            if rates[0][1]/(rates[0][0] + rates[0][1]) < rates[1][1]/(rates[1][0] + rates[1][1]):
                fout.write(f'{tree_path}\t{sig}\t{node_name}\t|\t{",".join(in_ref)}\t|\t{",".join(out_ref)}\n')
            else :
                fout.write(f'{tree_path}\t{sig}\t{node_name}\t|\t{",".join(out_ref)}\t|\t{",".join(in_ref)}\n')

if __name__ == "__main__":
    main()