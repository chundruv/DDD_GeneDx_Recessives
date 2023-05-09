import numpy as np

def is_inherited_proband(variant, child_id, ped, samples, gt_types):
    dad_id=ped[child_id][0]
    mum_id=ped[child_id][1]
    child_geno=gt_types[np.where(np.isin(samples,child_id))[0]]
    dad_geno=gt_types[np.where(np.isin(samples,dad_id))[0]]
    mum_geno=gt_types[np.where(np.isin(samples,mum_id))[0]]
    if child_geno==3 or mum_geno==3 or dad_geno==3:
        return "NA"
    elif child_geno==1:
        if dad_geno>0 and mum_geno>0:
            return "unknown"
        elif dad_geno>0 and mum_geno==0:
            return "dad"
        elif dad_geno==0 and mum_geno>0:
            return "mum"
        elif dad_geno==0 and mum_geno==0:
            return "denovo"
    elif child_geno==2:
        if dad_geno>0 and mum_geno>0:
            return "dad/mum"
        elif dad_geno==0 and mum_geno>0:
            return "dad/denovo"
        elif dad_geno>0 and mum_geno==0:
            return "mum/denovo"
        elif dad_geno==0 and mum_geno==0:
            return "denovo/denovo"
    elif child_geno==0:
        return "non_transmitted"

def is_inherited_parent(variant, parent_id, parents, ped ,samples, gt_types):
    child_geno=gt_types[np.where(np.isin(samples,parents[parent_id]))[0]]
    dad_geno=gt_types[np.where(np.isin(samples,ped[parents[parent_id]][0]))[0]]
    mum_geno=gt_types[np.where(np.isin(samples,ped[parents[parent_id]][1]))[0]]

    if child_geno==3:
        return "NA"
    elif child_geno>0:
        if child_geno==1 and dad_geno>0 and mum_geno>0:
            return "unclear"
        elif child_geno==1 and dad_geno>0 and mum_geno==0:
            return "dad"
        elif child_geno==1 and dad_geno==0 and mum_geno>0:
            return "mum"
        elif (child_geno==1 and (dad_geno==mum_geno==0)):
            return "denovo"
        elif child_geno==2 and dad_geno>0 and mum_geno>0:
            return "both"
        elif child_geno==2 and (dad_geno==mum_geno==0):
            return "double-denovo/likely error"
        else:
            return "inh/denovo"
    else:
        return "non_transmitted"

def get_inh(samples_pop_inds, parents, ped, samples, gt_types, variant, ind):
    pop=''
    child_inh='NA'
    parent_inh='NA'

    if ind in samples_pop_inds:
        pop=samples_pop_inds[ind]
        if ind in parents:
            parent_inh=is_inherited_parent(variant, ind, parents, ped, samples, gt_types)
            child_inh='NA'
        elif ind in ped:
            parent_inh='NA'
            child_inh=is_inherited_proband(variant, ind, ped, samples, gt_types)
    return (pop, child_inh, parent_inh)
