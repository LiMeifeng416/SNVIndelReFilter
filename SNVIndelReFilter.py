import argparse
import json
import re
#########################################################################################################################
def Judge_item(item, Select_paramer,head):
    S_value = {}
    index = head.split("\t")
    for i in range(0,len(index)):
        name=index[i]
        value = item.split("\t")[i]
        if "%" in value:
            value=float(value.strip("%"))/100
        elif value.isdigit():
            value=float(value)

        S_value[name]=value
    formats=Select_paramer.format_map(S_value)
    if eval(formats):
        result=True
    else:
        result=False
    return result


def out_file(dict,out_file,K_index):
        for key in dict:
            dter = dict[key][K_index].strip("%")
            if float(dter) > 100:
                dict[key][K_index] = "100.00%"
                dict[key][K_index - 1] = dict[key][K_index + 1]
            out_file.write("\t".join(dict[key]))

def Pass_judge(info,dict,K_index):
    idx = tuple(info[:5])
    if idx not in dict:
        dict[idx] = info
    else:
        dter = dict[idx][K_index].strip("%")
        dter2 = info[K_index].strip("%")
        if float(dter2) > float(dter):
            dict[idx] = info

def Dis_judge(info, Dic_union,dict, K_index):
    idx = tuple(info[:5])
    if (idx not in Dic_union) & (idx not in dict):
        dict[idx] = info
    elif (idx not in Dic_union) & (idx in dict):
        dter = dict[idx][K_index].strip("%")
        dter2 = info[K_index].strip("%")
        if float(dter2) > float(dter):
            dict[idx] = info
#########################################################################################################################
def SNVIndelReFilter(sample_flt_ann, sample_var_ann, sample_dis_ann, out_flt, out_dis1, out_dis2,json_file):
    glst=['Chr','Start','End','Ref','Alt','Tags','Freq','Amplicon']
    minus_set=["C-T","G-T"]
    plus_set=["G-A","C-A"]
    dict_dis_1={}
    dict_flt={}
    Dic_union={}

    with open(json_file) as fp:
        config=json.load(fp)
        Lable=config["Lable"]
        Tag_paramer=config["Tags"]
        Select_paramr=config["Select_set"]
        lab=config["Key"]

    with open(sample_flt_ann)as flt_ann, open(out_flt, 'w') as flt_2, open(sample_dis_ann) as dis_ann,open (out_dis1,'w') as dis_ann_1,open (out_dis2,'w') as dis_ann_2:
        H = next(flt_ann)
        next(dis_ann)
        K_index = H.split("\t").index(lab)
        index_tags = H.split("\t").index('Tags')
        flt_2.write(H)
        dis_ann_1.write(H)
        dis_ann_2.write(H)
        head = H.strip().split("\t")
        c, s, e, r, v, t, f, amp = [head.index(i) for i in glst]
        for temp in flt_ann:
            i = temp.split("\t")
            ref = i[r]
            var = i[v]
            mut = ref + "-" + var
            plus = int(i[-2])
            minus = int(i[-1])
            adp = int(i[9])
            freq = float(i[f].strip("%"))
            amplicon = i[amp].strip("*").split(";")

            if ((mut in minus_set and plus == 0) or (
                    mut in plus_set and minus == 0)) and adp <= 30 and freq <= 10 and len(amplicon) == 1:
                with open(sample_var_ann) as var_ann:
                    var_lst = ['Ref', 'Var', 'Total', 'SSBC', 'Amplicon']
                    head_var = next(var_ann)
                    head_var = head_var.strip().split("\t")
                    rv, vv, tv, sv, av = [head_var.index(m) for m in var_lst]
                    num = 0
                    for j in var_ann:
                        j = j.split("\t")
                        amplicon_var = j[av].strip("*")
                        if (amplicon_var == amplicon[0]) and (j[rv] == ref and j[vv] == var) and Lable=="Total":
                            freq_temp = float(j[tv].split("AF=")[1].split(";")[0]) * 100
                            ADP_temp = int(j[tv].split("AD=")[1].split(";")[0])
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2):
                                num += 1
                        elif (amplicon_var == amplicon[0]) and (j[rv] == ref and j[vv] == var) and Lable=="SSBC":
                            freq_temp = float(j[sv].split("AF=")[1].split(";")[0]) * 100
                            ADP_temp = int(j[sv].split("AD=")[1].split(";")[0])
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2):
                                num += 1

                if num > 2:
                    out_discrad = temp.replace("PASS", "Suspect")
                    tags = out_discrad.split("\t")[index_tags]
                    lable = True
                    for j in Tag_paramer:
                        if j in tags:
                            lable = False
                            dis_ann_2.write(out_discrad)
                            break
                    if lable:
                        Key = Judge_item(out_discrad, Select_paramr, H)
                        if Key:
                            info = out_discrad.split("\t")
                            Dis_judge(info,Dic_union,dict_dis_1,K_index)
                        else:
                            dis_ann_2.write(out_discrad)

                else:
                    Pass_judge(i,dict_flt,K_index)
            else:
                Pass_judge(i,dict_flt,K_index)

            Dic_union.update(dict_flt)

        for m in dis_ann:
            tags = m.split("\t")[index_tags]
            L = True
            for j in Tag_paramer:
                if j in tags:
                    L = False
                    dis_ann_2.write(m)
                    break
            if L:
                Key = Judge_item(m, Select_paramr, H)
                if Key:
                    info = m.split("\t")
                    idx = tuple(info[:5])
                    if (idx not in Dic_union) & (idx not in dict_dis_1):
                        dict_dis_1[idx] = info
                    elif (idx not in Dic_union) & (idx in dict_dis_1):
                        dter = dict_dis_1[idx][K_index].strip("%")
                        dter2 = info[K_index].strip("%")
                        if float(dter2) > float(dter):
                            dict_dis_1[idx] = info
                else:
                    dis_ann_2.write(m)

        out_file(dict_flt,flt_2,K_index)
        out_file(dict_dis_1,dis_ann_1, K_index)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SNVIndelReFilter")
    parser.add_argument("-i", "--flt", help="flt.ann file", required=True)
    parser.add_argument("-v", "--var", help="var.ann file", required=True)
    parser.add_argument("-d", "--dis", help="dis.ann", required=True)
    parser.add_argument("-of", "--outflt", help="out flt.ann2 file", required=True)
    parser.add_argument("-o1", "--simplify", help="simplify_suffix", required=True)
    parser.add_argument("-o2", "--other", help="other_suffix", required=True)
    parser.add_argument("-j", "--json", help="json file", required=True)
    args = parser.parse_args()
    SNVIndelReFilter(sample_flt_ann=args.flt, sample_var_ann=args.var, sample_dis_ann=args.dis, out_flt=args.outflt,
                     out_dis1=args.simplify, out_dis2=args.other,json_file=args.json)

