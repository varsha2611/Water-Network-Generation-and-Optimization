def plot_deviation_2(vals_of_replicas, vals_of_graph, metrics, figpath, jaccard_edges=None, title_infix='', seed=0, Gname=''):
#vals_of_graph could be a number (level 0) or a list (the same as the number of replicas)
    clean_names = {'num nodes': 'num nodes', 'num edges':'num edges', 'clustering':'clustering', 'average degree':'avg\ndegree', 'degree assortativity':'degree\nassortativity', 'degree connectivity':'degree\nconnectivity',
            'total deg*deg':'total deg*deg\nassortativity',
            's-metric':'s metric', 'mean ecc':'avg\neccentricity', 'num comps':'num comps', 'L eigenvalue sum':'L eigen-\nvalue sum',
            'average shortest path':'avg\ndistance', 'harmonic mean path':'harmonic avg\ndistance', 'avg flow closeness':'avg flow\ncloseness',
            'avg eigvec centrality':'avg eigenvec.\ncentrality', 'avg between. central.':'avg between.\ncentrality', 'modularity':'modularity'}

    multiple_models = type(vals_of_graph[0]) is list

    pylab.show()
    fig = pylab.figure()
    pylab.hold(True)
    num_of_metrics = len(metrics)
    med_vals = [np.median(vals_of_replicas[i]) for i in range(num_of_metrics)]
    avg_vals = [np.average(vals_of_replicas[i]) for i in range(num_of_metrics)]
    p25_vals = [np.percentile(vals_of_replicas[i],25) for i in range(num_of_metrics)]
    p75_vals = [np.percentile(vals_of_replicas[i],75) for i in range(num_of_metrics)]
    max_vals = [np.max(vals_of_replicas[i]) for i in range(num_of_metrics)]
    min_vals = [np.min(vals_of_replicas[i]) for i in range(num_of_metrics)]
    std_vals = [np.std(vals_of_replicas[i]) for i in range(num_of_metrics)]

    replica_stats = {'median_of_replicas':med_vals, 'avg_of_replicas':avg_vals, 'p25_of_replicas':p25_vals, 'p75_of_replicas':p75_vals, 'max_of_replicas':max_vals, 'min_of_replicas':min_vals, 'std_of_replicas':std_vals}

    normed_replica_vals = []
    avg_norms  = []
    print('Medians' + (' (average of model graphs)' if multiple_models else ''))
    print('-------')
    print('metric\t\tOriginalG\t\tReplicas')
    for met_num,metric in enumerate(metrics):
        try:
            model_val = np.average(vals_of_graph[met_num]) if multiple_models else vals_of_graph[met_num]
            print('%s\t\t%.5f\t\t%.5f'%(metric['name'],model_val,med_vals[met_num]))
        except:
            print('%\tserror'%metric['name'])
    for met_num,metric in enumerate(metrics):
        #handle error in original, 0 in original, error in one replica, error in all replicas
        nor_vals = []
        if multiple_models:
            assert len(vals_of_graph[met_num]) == len(vals_of_replicas[met_num])
            pruned_model_vals = [v for v in vals_of_graph[met_num] if v!=graphutils.METRIC_ERROR]
            if len(pruned_model_vals) > 0:
                v_graph = np.average(pruned_model_vals)
            else:
                v_graph = graphutils.METRIC_ERROR
        else:
            v_graph     = vals_of_graph[met_num]

        v_reps      = vals_of_replicas[met_num]
        if v_graph != graphutils.METRIC_ERROR:
            if v_graph != 0.0:
                nor_vals = [float(v)/v_graph for v in v_reps if v != graphutils.METRIC_ERROR]
            else:
                if v_reps != [] and np.abs(v_reps).sum() == 0.:
                    nor_vals.append(len(v_reps)*[1.0])
            pylab.plot(1.0, met_num, 'o', color='k', linewidth=2., label=Gname)
            pylab.text(x=.0, y=(met_num-2./len(metrics)), s='%.2e'%v_graph)
            #if type(v_graph) is int:
            #    pylab.text(x=.0, y=(met_num-2./len(metrics)), s=str(v_graph))
            #else:
            #    pylab.text(x=.0, y=(met_num-2./len(metrics)), s='%.3f'%v_graph)
            nor_vals = np.array(nor_vals)
            normed_replica_vals.append(nor_vals)
            if len(nor_vals) >0:
                bp = pylab.boxplot(nor_vals, positions=[met_num], vert=1, widths=0.5, patch_artist=True)
                for patch in (bp['boxes']):
                    patch.set_facecolor('lightblue')
                ax.yaxis.grid(True)
                ax.set_xticks([y + 1 for y in range(len(all_data))], )
                ax.set_xlabel('xlabel')
                ax.set_ylabel('ylabel')
                plt.setp(axes, xticks=[y + 1 for y in range(len(nor_vals))])
                plt.show()
            if (nor_vals == graphutils.METRIC_ERROR).any():
                    val_str= r'undefined'
                    avg_norm = -np.inf
                elif np.abs(nor_vals).sum() < 1000:
                    avg_norm = np.average(nor_vals)
                    val_str= r'$%.2f$'%np.average(nor_vals) if latex_available else r'%.2f'%avg_norm
                else:
                    avg_norm = np.inf
                    val_str= r'$\gg0$'if latex_available else r'>>0'
                avg_norms.append(avg_norm)
            else:
                val_str = r'undefined'
                avg_norms.append(None)
        else:
            val_str = r'undefined'
            normed_replica_vals.append([None, None])
            avg_norms.append(None)
        pylab.text(x=1.74, y=(met_num-2./len(metrics)), s=val_str)
    try:
        pylab.yticks(list(range(num_of_metrics)), [clean_names.get(met['name'], met['name']) for met in metrics], rotation=0)
        if multiple_models:
            pylab.xlabel(r'Relative to mean of coarse networks', rotation=0, fontsize='20')#, x=0.1)
        else:
            pylab.xlabel(r'Relative to real network', rotation=0, fontsize='20')#, x=0.1)
        #pylab.title(G.name)
        #pylab.legend(loc='best')
        max_axis = 2
        pylab.xlim(-0.02,max_axis)
        pylab.ylim(-1.0,len(metrics))
        pylab.text(x=0.00, y=len(metrics)+0.05, s='Template\ngraph', va='bottom')
        pylab.text(x=1.650, y=-1.05, s='Median of\nreplicas', va='top')
        if jaccard_edges != None:
            pylab.text(x=0.30, y=len(metrics)+0.05, s='(Jaccard=%.3f)'%jaccard_edges, va='bottom')
            #pylab.text(x=-0.30, y=len(metrics)*(-0.15), s='E[EdgeJaccard]=%.3f'%jaccard_edges, ha='right', va='top')

        fig.subplots_adjust(left=0.17, right=0.95)

    except Exception as inst:
        print('Warning: could not save stats figure '+figpath + ':\n'+str(inst))
        exc_traceback = sys.exc_info()[2]
        print(str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n'))

    replica_stats['normed_replica_vals'] = normed_replica_vals
    replica_stats['avg_norm_of_replicas'] = avg_norms

    mean_rel_errors = []
    mean_relstd_errors = []
    for met_i in range(num_of_metrics):
        normed_vals = normed_replica_vals[met_i]
        if graphutils.METRIC_ERROR in normed_vals or len(normed_vals) == 1:
            mean_rel_errors.append(None)
            mean_relstd_errors.append(None)
            continue
        rel_error_ar = [v - 1.0 for v in normed_vals if v != None]
        if len(rel_error_ar) == 0:
            rel_error_ar = [graphutils.METRIC_ERROR,  graphutils.METRIC_ERROR]
        mean_rel_errors.append(np.average(rel_error_ar))
        mean_relstd_errors.append(np.average(rel_error_ar)/(1E-20 + np.std(rel_error_ar)))

    replica_stats['mean_rel_errors'] = mean_rel_errors
    replica_stats['mean_relstd_errors'] = mean_relstd_errors
    try:
        replica_stats['mean_mean_error']    = np.average(mean_rel_errors)     #the grand stat
        replica_stats['mean_mean_errorstd'] = np.average(mean_relstd_errors)  #the grand stat
    except:
        replica_stats['mean_mean_error']    = None
        replica_stats['mean_mean_errorstd'] = None

    return replica_stats, figpath
