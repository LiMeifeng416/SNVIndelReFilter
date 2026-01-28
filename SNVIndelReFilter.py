#-*-coding:UTF-8-*-
import argparse
"""
Name:   SNVIndelReFilter
Version:    0.1.0
"""
def SNVIndelReFilter(sample_flt_ann,sample_var_ann,sample_dis_ann,out_flt_ann,out_dis_ann,lable):
    glst=['Chr','Start','End','Ref','Alt','Tags','Freq','Amplicon']
    minus_set=["C-T","G-T"]
    plus_set=["G-A","C-A"]
    with open(sample_flt_ann)as flt_ann, open(out_flt_ann, 'w') as flt_2, open(sample_dis_ann) as dis_ann,open (out_dis_ann,'w') as dis_ann_2:
        for m in dis_ann:
            dis_ann_2.write(m)
        head = next(flt_ann)
        flt_2.write(head)
        head = head.strip().split("\t")
        c, s, e, r, v, t, f, amp = [head.index(i) for i in glst]
        for temp in flt_ann:
            i = temp.strip().split("\t")
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
                        if (amplicon_var == amplicon[0]) and (j[rv] == ref and j[vv] == var) and lable=="Total":
                            freq_temp = float(j[tv].split("AF=")[1].split(";")[0]) * 100
                            ADP_temp = int(j[tv].split("AD=")[1].split(";")[0])
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2):
                                num += 1
                        elif (amplicon_var == amplicon[0]) and (j[rv] == ref and j[vv] == var) and lable=="SSBC":
                            freq_temp = float(j[sv].split("AF=")[1].split(";")[0]) * 100
                            ADP_temp = int(j[sv].split("AD=")[1].split(";")[0])
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2):
                                num += 1
                        else:
                            next
                if num > 2:
                    out_discrad = temp.replace("PASS", "Suspect")
                    dis_ann_2.write(out_discrad)
                else:
                    flt_2.write(temp)
            else:
                flt_2.write(temp)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SNVIndelReFilter")
    parser.add_argument("-i", "--flt",help="flt.ann file",required=True)
    parser.add_argument("-v", "--var", help="var.ann file",required=True)
    parser.add_argument("-d", "--dis", help="dis.ann",required=True)
    parser.add_argument("-of", "--outflt", help="out flt.ann2 file", required=True)
    parser.add_argument("-od", "--outdis", help="out dis.ann2 file", required=True)
    parser.add_argument("-k", "--key", help="SSBC or Total", required=True)
    args = parser.parse_args()
    SNVIndelReFilter(sample_flt_ann=args.flt,sample_var_ann=args.var,sample_dis_ann=args.dis,out_flt_ann=args.outflt,out_dis_ann=args.outdis,lable=args.key)
