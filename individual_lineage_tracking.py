import click, ete3, pandas as pd, numpy as np, collections, json, sys
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable  # Added for CMH test

#before_FMT & non-responder are counted as disease
#healthy & responder are counted as healthy
conv = {'before_FMT':'D', 'non-responder':'D',
        'healthy':'H', 'responder':'H'}

h_codes = {'H': 0, 'D': 1}
c_codes = {'rCDI':0, 'IBS':1, 'LUAD':2, 'MEL': 3}
cohort_names = {v: k for k, v in c_codes.items()}  # Reverse mapping for cohorts

def get_optimal_cut(tre, samples) :
    tips = {}
    for tip in tre.get_leaves() :
        info = tip.name.split('|')
        if info[0] in samples :
            tips[tip.name] = samples[info[0]]
    
    branches = []
    for n in tre.traverse('postorder') :
        if n.is_leaf() :
            if n.name in tips :
                n.d = set([n.name])
            else :
                n.d = set([])
        else :
            n.d = set([d for c in n.children for d in c.d])
            if len(n.d) > 1 :
                ingroup = collections.defaultdict(list)
                outgroup = collections.defaultdict(list)
                for t, info in tips.items() :
                    if t in n.d :
                        ingroup[info].append(t.split('|')[0])
                    else :
                        outgroup[info].append(t.split('|')[0])
                if len(ingroup) > 1 and len(outgroup) > 1 :
                    # Disease counts (for output only)
                    d_h = [
                        np.bincount([h_codes[x[2]] for x in ingroup.keys()], minlength=len(h_codes)),
                        np.bincount([h_codes[x[2]] for x in outgroup.keys()], minlength=len(h_codes))
                    ]
                    if np.sum(d_h[0]) < 0.1*np.sum(d_h) or np.sum(d_h[0]) > 0.9*np.sum(d_h):
                        continue
                    # Cohort counts
                    d_c = [
                        np.bincount([c_codes[x[0]] for x in ingroup.keys()], minlength=len(c_codes)),
                        np.bincount([c_codes[x[0]] for x in outgroup.keys()], minlength=len(c_codes))
                    ]
                    if np.sum(np.min(d_c, 0) >= 0.1*np.sum(d_c, 0)) < 0.5*len(d_c[0]):
                        continue

                    # Handle cohort test (original Fisher's exact for cohort distribution)
                    res_h = fisher_exact(d_h)
                    health_fisher = res_h[1] if isinstance(res_h, tuple) else res_h
                    res_c = fisher_exact(d_c)
                    cohort_fisher = res_c[1] if isinstance(res_c, tuple) else res_c
                    
                    # ===== NEW: CMH TEST FOR DISEASE (ADJUSTED FOR COHORT) =====
                    # Build stratified tables by cohort
                    cmh_tables = []
                    for cohort_idx in range(len(c_codes)):
                        cohort_name = cohort_names[cohort_idx]
                        in_H = 0
                        in_D = 0
                        out_H = 0
                        out_D = 0
                        
                        # Count disease status in ingroup for this cohort
                        for info in ingroup.keys():
                            if info[0] == cohort_name:
                                if info[2] == 'H': in_H += 1
                                else: in_D += 1
                                
                        # Count disease status in outgroup for this cohort
                        for info in outgroup.keys():
                            if info[0] == cohort_name:
                                if info[2] == 'H': out_H += 1
                                else: out_D += 1
                        
                        # Only include cohorts with data
                        if in_H + in_D + out_H + out_D > 0:
                            cmh_tables.append([[in_H, in_D], [out_H, out_D]])
                    
                    # Run CMH test if valid tables exist
                    cmh_pvalue = 1.0  # Default if test fails
                    if len(cmh_tables) > 0:
                        try:
                            st = StratifiedTable(cmh_tables, shift_zeros=True)
                            cmh_result = st.test_null_odds()
                            cmh_pvalue = cmh_result.pvalue
                        except Exception as e:
                            print(f"CMH error: {e}", file=sys.stderr)
                    # ===== END CMH TEST =====
                    
                    
                    branches.append([
                        cmh_pvalue,            # Disease p-value (CMH-adjusted)
                        health_fisher,
                        cohort_fisher,         # Cohort p-value (Fisher's)
                        n.name, 
                        np.array(d_h).tolist(),
                        np.array(d_c).tolist(), 
                        list(set([ss for s in ingroup.values() for ss in s])), 
                        list(set([ss for s in outgroup.values() for ss in s]))
                    ])

    return min(branches) if branches else None

@click.command()
@click.option('-m', '--metadata')
@click.option('-n', '--nwk')
def main(metadata, nwk) :
    meta = pd.read_csv(metadata, sep='\t', header=0)
    samples = {}
    for id, cohort, individual, day, dtype in meta[['ID', 'Cohort', 'individual', 'Day', 'Disease type']].values :
        samples[id] = (cohort, individual, conv[dtype])

    tre = ete3.Tree(nwk, format=1)
    data = get_optimal_cut(tre, samples)
    if data:
        print(f'{nwk}\t{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\t{json.dumps(data[4]).replace(" ", "")}\t{json.dumps(data[5]).replace(" ", "")}\t{json.dumps(data[6]).replace(" ", "")}\t{json.dumps(data[7]).replace(" ", "")}')

if __name__ == '__main__' :
    main()
