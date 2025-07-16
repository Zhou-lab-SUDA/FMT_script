import click, ete3, pandas as pd, numpy as np, collections


def get_distance(tre, samples) :
    individual_pairs = collections.defaultdict(lambda : collections.defaultdict(dict))
    for n in tre.traverse('postorder') :
        if n.is_leaf() :
            info = n.name.split('|')
            if info[0] in samples :
                n.d = [[n.name, info[0], n.dist, 0.]]
            else :
                n.d = []
        else :
            for i, c1 in enumerate(n.children) :
                for j, c2 in enumerate(n.children[:i]) :
                    if len(c1.d) and len(c2.d) :
                        for t1 in c1.d :
                            for t2 in c2.d :
                                if (t1[1] != t2[1]) and (samples[t1[1]][0] == samples[t2[1]][0]) and samples[t1[1]][1] != samples[t2[1]][1] :
                                    diff_day = tuple(sorted([samples[t1[1]][1], samples[t2[1]][1]]))
                                    dist = t1[3] + t2[3] + (t1[2] if 'low_qual' in t1[0] else t1[2]/4) + (t2[2] if 'low_qual' in t2[0] else t2[2]/4)
                                    if diff_day not in individual_pairs[samples[t1[1]][0][0]][samples[t1[1]][0][1]] or individual_pairs[samples[t1[1]][0][0]][samples[t1[1]][0][1]][diff_day] > dist :
                                        individual_pairs[samples[t1[1]][0][0]][samples[t1[1]][0][1]][diff_day] = dist
            n.d = [[d[0], d[1], d[2], d[3] + n.dist] for c in n.children for d in c.d]
    return individual_pairs

@click.command()
@click.option('-m', '--metadata')
@click.option('-n', '--nwk')
@click.option('-p', '--prefix')
def main(prefix, metadata, nwk) :
    meta = pd.read_csv(metadata, sep='\t', header=0)
    samples = {}
    for id, cohort, individual, day, dtype in meta.loc[meta['Disease type'] != 'before_FMT', ['ID', 'Cohort', 'individual', 'Day', 'Disease type']].values :
        samples[id] = [(cohort, individual), day, dtype, cohort]

    tre = ete3.Tree(nwk, format=1)
    data = get_distance(tre, samples)
    print(f'Prefix,Cohort,Num_individuals,delta_Date,mean_Persistence,median,2.5%,25%,75%,97.5%')
    for cohort, individuals in sorted(data.items()) :
        individuals = np.array(list(individuals.values()))
        indicies = np.random.choice(len(individuals), len(individuals)*3000).reshape([3000, len(individuals)])
        results = collections.defaultdict(list)
        for idx in indicies :
            res = collections.defaultdict(list)
            for idv in individuals[idx] : 
                res[999999999].append(min(idv.values())<=0.001)
                for dates, dist in idv.items() :
                    dc = np.power(2, int(np.log2(dates[1] - dates[0])))
                    res[dc].append(dist <= 0.001)
            for dates, dists in res.items() :
                results[dates].append(np.average(dists))
        for dates, vals in sorted(results.items()) :
            qs = np.quantile(vals, [0.025, 0.25, 0.5, 0.75, 0.975])
            if dates == 999999999 :
                dates = 'SUM'
            print(f'{prefix},{cohort},{len(individuals)},{dates},{np.mean(vals):.2f},{qs[2]:.2f},{qs[0]:.2f},{qs[1]:.2f},{qs[3]:.2f},{qs[4]:.2f}')



if __name__ == '__main__' :
    main()