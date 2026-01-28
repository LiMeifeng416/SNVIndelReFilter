import argparse
import json
import re
from sys import argv
from smart_open import open
"""

History:
Version 2.0.0 
Version 2.1.0 Data:2023.12.14
1.增加对dis.ann文件中的损伤位点进行判断
2.将多次出现的位点标记符号由“*”变为“#”
3.将判断同批次中3个及以上相同非SNP位点的范围由PASS列表中的位点调整至PASS&discard中的所有位点
#未用：4.增加对整批样本fusion结果的注释，同批次中多次检出的融合添加“#”标记。

"""
#########################################################################################################################
##计算筛选select_set表达式是否满足条件
def Judge_item(item, Select_paramer,head):
    S_value = {}
    index = head.strip("\n").split("\t")
    for i in range(0,len(index)):
        name=index[i]
        value = item[i]
        if "Comment" not in name:
            if "%" in value:
                value=float(value.strip("%"))/100
            elif value.isdigit():
                value=float(value)
            elif "'" in value:
                value=value.replace("'","")
            S_value[name]=value
    formats=Select_paramer.format_map(S_value)
    if eval(formats):
        result=True
    else:
        result=False
    return result

##将储存在dict中的突变items，写出至对应的文件是否存在突变对应的频率超出100%（bug），若存在则频率按照100%输出。
def out_file(dict,out_file,K_index):
        for key in dict:
            dter = dict[key][K_index].strip("%")
            if float(dter) > 100:
                dict[key][K_index] = "100.00%"
                dict[key][K_index - 1] = dict[key][K_index + 1]
            out_file.write("\t".join(dict[key]))

##判断PASS_dict中是否存在相同的突变，若存在则将频率更高的items输出
def Pass_judge(info,dict,K_index):
    idx = tuple(info[:5])
    if idx not in dict:
        dict[idx] = info
    else:
        dter = dict[idx][K_index].strip("%")
        dter2 = info[K_index].strip("%")
        if float(dter2) > float(dter):
            dict[idx] = info

##判断DIS_dict中是否存在相同的突变（并且不在PASS_dict中），若存在则将频率更高的items输出
def Dis_judge(info, Dic_union,dict, K_index):
    idx = tuple(info[:5])
    if (idx not in Dic_union) & (idx not in dict):
        dict[idx] = info
    elif (idx not in Dic_union) & (idx in dict):
        dter = dict[idx][K_index].strip("%")
        dter2 = info[K_index].strip("%")
        if float(dter2) > float(dter):
            dict[idx] = info


##判断位点是否为损伤位点
###判断PASS列表位点损失情况
def Judge_damage(var_ann_file,flt_ann_file,dict_flt,dict_dis,K_index,head,Lable):
    glst=['Chr','Start','End','Ref','Alt','Tags','Freq','Amplicon','AltDepth','Plus','Minus']
    minus_set=["C-T","G-T"]
    plus_set=["G-A","C-A"]
    with open(flt_ann_file) as flt_ann:
        next(flt_ann)
        c, s, e, r, v, t, f, amp,ad,p,m = [head.index(i) for i in glst]
        for temp in flt_ann:
            i = temp.split("\t")
            ref = i[r]
            var = i[v]
            mut = ref + "-" + var
            plus = int(i[p])
            minus = int(i[m])
            adp = int(i[ad])
            freq = float(i[f].strip("%"))
            amplicon = i[amp].strip("*").split(";")
            ##判断突变是否为损伤类型的突变
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
                    info = out_discrad.split("\t")
                    Dis_judge(info,dict_flt,dict_dis,K_index)
                else:
                    Pass_judge(i,dict_flt,K_index)
            else:
                Pass_judge(i,dict_flt,K_index)

###判断dis列表位点损失情况，同时判断dis中是否存在相同突变检出（同时PASS中无检出）
def Judge_damage_dis(var_ann_file,dis_ann_file,dict_flt,dict_dis,K_index,head,Lable):
    glst=['Chr','Start','End','Ref','Alt','Tags','Freq','Amplicon','AltDepth','Plus','Minus']
    minus_set=["C-T","G-T"]
    plus_set=["G-A","C-A"]
    with open(dis_ann_file) as dis_ann:
        next(dis_ann)
        c, s, e, r, v, t, f, amp,ad,p,m = [head.index(i) for i in glst]
        for temp in dis_ann:
            i = temp.split("\t")
            ref = i[r]
            var = i[v]
            mut = ref + "-" + var
            plus = int(i[p])
            minus = int(i[m])
            adp = int(i[ad])
            freq = float(i[f].strip("%"))
            amplicon = i[amp].strip("*").split(";")
            ##判断突变是否为损伤类型的突变
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
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2) :
                                num += 1
                        elif (amplicon_var == amplicon[0]) and (j[rv] == ref and j[vv] == var) and Lable=="SSBC":
                            freq_temp = float(j[sv].split("AF=")[1].split(";")[0]) * 100
                            ADP_temp = int(j[sv].split("AD=")[1].split(";")[0])
                            if ADP_temp >= (adp / 2) or 1 <= freq_temp <= 3 or freq_temp>=(freq/2):
                                num += 1
                if num >= 2:
                    Tags_new="Suspect;"+i[t]
                    #out_discrad = temp.replace("PASS", "Suspect")
                    #info = out_discrad.split("\t")
                    i[t]=Tags_new
                    Dis_judge(i,dict_flt,dict_dis,K_index)
                else:
                    Dis_judge(i,dict_flt,dict_dis,K_index)
            else:
                Dis_judge(i,dict_flt,dict_dis,K_index)



##判断tags标签是否满足要求
def Judge_tags(item_tags,Tag_paramer):
    item_tags_list=item_tags.split(";")
    L = True
    for j in Tag_paramer:
        for item_tag in item_tags_list:
            if j == item_tag:
                L = False
                break
    return L

#########################################################################################################################
def SNVIndelReFilter_Sample(sample_flt_ann, sample_var_ann, sample_dis_ann, out_flt, out_dis1, out_dis2,json_file):
    glst=['Chr','Start','End','Ref','Alt','Tags','Freq','Amplicon']
    minus_set=["C-T","G-T"]
    plus_set=["G-A","C-A"]
    dict_flt={}
    dict_dis={}
    dict_dis_1={}
    dict_dis_2={}
    Low_lable=['Low_Freq','Low_Dep','Low_ADP']
    with open (json_file) as fp, open(out_flt, 'w') as flt_2, open(sample_dis_ann) as dis_ann,open (out_dis1,'w') as dis_ann_1,open (out_dis2,'w') as dis_ann_2:
        
        config=json.load(fp)
        Lable=config["Lable"]
        Tag_paramer_1=config["Tags"]
        Select_paramr_1=config["Select_set"]
        lab=config["Key"]
        Tag_paramer_2=config["QTags"]
        Select_paramr_2=config["QSelect_set"]


        H = next(dis_ann)
        K_index = H.split("\t").index(lab)
        index_tags = H.split("\t").index('Tags')
        flt_2.write(H)
        dis_ann_1.write(H)
        dis_ann_2.write(H)
        head = H.strip().split("\t")

        Judge_damage(sample_var_ann,sample_flt_ann,dict_flt,dict_dis,K_index,head,Lable)
        Judge_damage_dis(sample_var_ann,sample_dis_ann,dict_flt,dict_dis,K_index,head,Lable)

        for key in dict_dis.keys():
            Tags=dict_dis[key][index_tags]
            Value=dict_dis[key]
            Tags_list=Tags.split(";")
            if set(Tags_list).isdisjoint(Tag_paramer_2):
                lable_1=Judge_tags(Tags,Tag_paramer_1)
                lable_2=Judge_item(Value, Select_paramr_1,H)
                if lable_1 and lable_2:
                    dict_dis_1[key]=Value
                else:
                    dict_dis_2[key]=Value
            else:
                if set(Tags_list).issubset(Tag_paramer_2):
                    lable_1=Judge_item(Value, Select_paramr_2,H)
                    if lable_1:
                        dict_dis_1[key]=Value
                    else:
                        dict_dis_2[key]=Value
                else:
                    lable_1=Judge_tags(Tags,Tag_paramer_1)
                    lable_2=Judge_item(Value, Select_paramr_2,H)
                    lable_3=Judge_tags(Tags,Low_lable)
                    if lable_1 and lable_2 and lable_3:
                        dict_dis_1[key]=Value
                    else:
                        dict_dis_2[key]=Value
            

        out_file(dict_flt,flt_2,K_index)
        out_file(dict_dis_1,dis_ann_1, K_index)
        out_file(dict_dis_2,dis_ann_2, K_index)

############################################################################################################################################
##获取在指定文件范围内，每个位点被检出的次数
def read(file,dic):
    with open(file) as file_ann:
        head=next(file_ann)
        index=head.split("\t").index("CDSChange")
        Tag_index=head.split("\t").index("Tags")
        for i in file_ann:
            CDS=i.split("\t")[index]
            tag=i.split("\t")[Tag_index]
            if CDS not in dic or "Polymorphism" in tag :
                dic[CDS]=1
            else:
                dic[CDS]+=1

###获取每个样本在PASS中检出的热点突变的次数
def read_Hotspot(sample_name,file,dic,hot_type):
    dic[sample_name]=0
    with open(file) as file_ann:
        head=next(file_ann)
        index=head.split("\t").index("Tags")
        for i in file_ann:
            temp=i.split("\t")
            if Judge_item(temp,hot_type,head):
                dic[sample_name]+=1

#####对dis列表中的位点进行重过滤，在同批次中多次检出的位点加上标记。
def dis_Rewrite(file,dic,dis_out_file,num):
    with open(file) as file_ann,open(dis_out_file,'w') as dis_out_ann:
        head=next(file_ann)
        CDS_index=head.split("\t").index("CDSChange")
        Tag_index=head.split("\t").index("Tags")
        dis_out_ann.write(head)
        for i in file_ann:
            if dic[i.split("\t")[CDS_index]]<num:
                dis_out_ann.write(i)
            else:
                tag=i.split("\t")[Tag_index]
                replace_tag="!"+tag
                i=i.replace(tag,replace_tag)
                dis_out_ann.write(i)
            
###Pass列表中存在同一个突变多次检出，或者同一个样本检出多个HotSpot的位点被标记#


def flt_Rewrite(flt_file,dic,dic_HotSpot,flt_out_file,sample,num_p,num_h,hot_type):
    with open(flt_file) as flt_ann,open(flt_out_file,'w') as flt_out_ann:
        head=next(flt_ann)
        CDS_index=head.split("\t").index("CDSChange")
        Tag_index=head.split("\t").index("Tags")
        flt_out_ann.write(head)
        for i in flt_ann:
            temp=i.split("\t")
            if Judge_item(temp,hot_type,head):
                if dic[i.split("\t")[CDS_index]]<num_p and dic_HotSpot[sample]<num_h:
                    flt_out_ann.write(i)
                elif dic[i.split("\t")[CDS_index]]<num_p and dic_HotSpot[sample]>=num_h:
                    tag=i.split("\t")[Tag_index]
                    replece_tag="#"+tag
                    i=i.replace(tag,replece_tag)
                    flt_out_ann.write(i)
                elif dic[i.split("\t")[CDS_index]]>=num_p and dic_HotSpot[sample]>=num_h:
                    tag=i.split("\t")[Tag_index]
                    replece_tag="#!"+tag
                    i=i.replace(tag,replece_tag)
                    flt_out_ann.write(i)
                elif dic[i.split("\t")[CDS_index]]>=num_p and dic_HotSpot[sample]<num_h:
                    tag=i.split("\t")[Tag_index]
                    replece_tag="!"+tag
                    i=i.replace(tag,replece_tag)
                    flt_out_ann.write(i)
            else:
                if dic[i.split("\t")[CDS_index]]<num_p:
                    flt_out_ann.write(i)
                else:
                    tag=i.split("\t")[Tag_index]
                    replece_tag="!"+tag
                    i=i.replace(tag,replece_tag)
                    flt_out_ann.write(i)
#########################
###对fusion文件进行判断，同批次样本中检出多个相同融合型，且其中存在强阳融合，其余拷贝数相差10倍及以上的融合标记#

###确定所有检出融合中拷贝数最大的情况
def read_fus(fus_file,dic):
	with open(fus_file) as fus_ann:
		head=next(fus_ann)
		copy_index=head.split("\t").index("Copies")
		fusion_index=head.split("\t").index("Fusion")
		for i in fus_ann:
			fusion=i.split("\t")[fusion_index]
			copies=i.split("\t")[copy_index]
			if fusion not in dic:
				dic[fusion]=copies
			else:
				if int(copies) > int(dic[fusion]):
					dic[fusion]=copies

####对拷贝数相差10倍及以上的融合进行注释“#”
def Rewrite_fus(fus_file,dic,out_fus_file):
	with open(fus_file) as fus_ann,open(out_fus_file,'w') as out_fus_ann:
		head=next(fus_ann)
		copy_index=head.split("\t").index("Copies")
		fusion_index=head.split("\t").index("Fusion")
		out_fus_ann.write(head)
		for i in fus_ann:
			fusion=i.split("\t")[fusion_index]
			copies=int(i.split("\t")[copy_index])
			if copies <= int(dic[fusion])/10:
				replace_fusion="#"+fusion
				i=i.replace(fusion,replace_fusion)
				out_fus_ann.write(i)
			else:
				out_fus_ann.write(i)

##########################################################################


def SNVIndelReFilter_Summary(filelist_file,flt_suffix,dis_suffix,discard_suffix,out_flt_suffix,out_dis_suffix,dirs,json_file):
	Dic={}
	Dic_flt={}
	Dic_HotSpot={}
	control_Dic={}
	control_Dic_flt={}
	control_Dic_HotSpot={}
	Dic_fus={}

	with open (json_file) as fp:
		config=json.load(fp)
		num_dis=config["number_dis"]
		num_pass=config["number_pass"]
		num_hot=config["number_HotSpot"]
		hot_type=config["HotSpot_type"]
		if 'fusion_suffix' in config:
			fusion_suffix=config["fusion_suffix"]

			with open(filelist_file) as FILELIST:
				for i in FILELIST:
					sample=i.split(" ")[0]
					flt_fus=dirs+sample+"/"+sample+".flt"+fusion_suffix
					dis_fus=dirs+sample+"/"+sample+".dis"+fusion_suffix
					read_fus(flt_fus,Dic_fus)
					read_fus(dis_fus,Dic_fus)
			with open(filelist_file) as FILELIST:
				for i in FILELIST:
					sample=i.split(" ")[0]
					flt_fus=dirs+sample+"/"+sample+".flt"+fusion_suffix
					out_flt_fus=dirs+sample+"/"+sample+".flt"+fusion_suffix+"2"
					dis_fus=dirs+sample+"/"+sample+".dis"+fusion_suffix
					out_dis_fus=dirs+sample+"/"+sample+".dis"+fusion_suffix+"2"
					Rewrite_fus(flt_fus,Dic_fus,out_flt_fus)
					Rewrite_fus(dis_fus,Dic_fus,out_dis_fus)

	####获取不同位点的重复次数
	with open(filelist_file) as FILELIST:
		for i in FILELIST:
			length=len(i.split(" "))
			if length==3:
				sample=i.split(" ")[0]
				flt_file=dirs+sample+"/"+sample+flt_suffix
				dis_file=dirs+sample+"/"+sample+dis_suffix
				if discard_suffix:
					discard_file=dirs+sample+"/"+sample+discard_suffix
					read(discard_file,Dic)
				read(flt_file,Dic)
				read(flt_file,Dic_flt)
				read(dis_file,Dic)
				read_Hotspot(sample,flt_file,Dic_HotSpot,hot_type)
			else:
				sample=i.split(" ")[0]
				control_sample=i.split(" ")[3]
				flt_file=dirs+sample+"/"+sample+flt_suffix
				dis_file=dirs+sample+"/"+sample+dis_suffix
				if discard_suffix:
					discard_file=dirs+sample+"/"+sample+discard_suffix
					read(discard_file,Dic)
				read(flt_file,Dic)
				read(flt_file,Dic_flt)
				read(dis_file,Dic)
				read_Hotspot(sample,flt_file,Dic_HotSpot,hot_type)
				
				control_flt_file=dirs+sample+"/"+control_sample+flt_suffix
				control_dis_file=dirs+sample+"/"+control_sample+dis_suffix
				if discard_suffix:
					control_discard_file=dirs+sample+"/"+control_sample+discard_suffix
					read(control_discard_file,control_Dic)
				read(control_flt_file,control_Dic)
				read(control_flt_file,control_Dic_flt)
				read(control_dis_file,control_Dic)
				read_Hotspot(sample,control_flt_file,control_Dic_HotSpot,hot_type)
	####根据获取到的不同位点的重复次数，对文件进行重新读写
	with open(filelist_file) as FILELIST:
		for i in FILELIST:
			length=len( i.split(" ") )
			if length==3:
				sample=i.split(" ")[0]
				flt_file=dirs+sample+"/"+sample+flt_suffix
				dis_file=dirs+sample+"/"+sample+dis_suffix
				dis_file_out=dirs+sample+"/"+sample+out_dis_suffix
				flt_file_out=dirs+sample+"/"+sample+out_flt_suffix
				dis_Rewrite(dis_file,Dic,dis_file_out,num_dis)
				flt_Rewrite(flt_file,Dic,Dic_HotSpot,flt_file_out,sample,num_pass,num_hot,hot_type)
			else:
				sample=i.split(" ")[0]
				flt_file=dirs+sample+"/"+sample+flt_suffix
				dis_file=dirs+sample+"/"+sample+dis_suffix
				dis_file_out=dirs+sample+"/"+sample+out_dis_suffix
				flt_file_out=dirs+sample+"/"+sample+out_flt_suffix
				dis_Rewrite(dis_file,Dic,dis_file_out,num_dis)
				flt_Rewrite(flt_file,Dic,Dic_HotSpot,flt_file_out,sample,num_pass,num_hot,hot_type)

				control_sample=i.split(" ")[3]
				control_flt_file=dirs+sample+"/"+control_sample+flt_suffix
				control_dis_file=dirs+sample+"/"+control_sample+dis_suffix
				control_dis_file_out=dirs+sample+"/"+control_sample+out_dis_suffix
				control_flt_file_out=dirs+sample+"/"+control_sample+out_flt_suffix
				dis_Rewrite(control_dis_file,control_Dic,control_dis_file_out,num_dis)
				flt_Rewrite(control_flt_file,control_Dic,control_Dic_HotSpot,control_flt_file_out,sample,num_pass,num_hot,hot_type)


###########################################################################################################################################
def argparser():
    desc='''Samples and Summary SNVIndelReFilter '''
    prog='SNVIndelReFilter'
    version='v2.1.0'
    parser=argparse.ArgumentParser(description=desc,prog=prog,add_help=False)
    required=parser.add_argument_group('optional')
    subparsers=parser.add_subparsers(help='sub-command help')
    j = subparsers.add_parser('Sample', help='ReFilter for samples')
    j.add_argument("-i", "--flt", help="flt.ann file", required=True)
    j.add_argument("-v", "--var", help="var.ann file", required=True)
    j.add_argument("-d", "--dis", help="dis.ann", required=True)
    j.add_argument("-of", "--outflt", help="out flt file", required=True)
    j.add_argument("-o1", "--simplify", help="simplify discard file", required=True)
    j.add_argument("-o2", "--other", help="check discard file", required=True)
    j.add_argument("-j", "--json", help="json file", required=True)

    f= subparsers.add_parser('Summary', help='ReFilter for Summary')
    f.add_argument("-i", "--filelist", help="FILELIST file", required=True)
    f.add_argument("-f", "--flt", help="suffix of flt file", required=True)
    f.add_argument("-d", "--dis", help="suffix of simplify discard file", required=True)
    f.add_argument("-dis", "--discard", help="suffix of check discard file", default=False)
    f.add_argument("-of", "--outflt", help="suffix of out flt file", required=True)
    f.add_argument("-od", "--outdis", help="suffix of out simplipy discard file", required=True)
    f.add_argument("-w", "--dir", help="dir", required=True)
    f.add_argument("-j", "--json", help="json file", required=True)

    help_group = parser.add_argument_group('Help optional')
    help_group.add_argument('-v', '--version', action='version', version="%s_%s" %(prog, version))
    help_group.add_argument('-h', '--help', action='help',help='show this help message and exit')
    return parser.parse_args()
    
if __name__ == "__main__":
    args = argparser()
    tp =argv[1]
    if tp =='Sample':
        sample_flt_ann=args.flt
        sample_var_ann=args.var
        sample_dis_ann=args.dis
        out_flt=args.outflt
        out_dis1=args.simplify
        out_dis2=args.other
        json_file=args.json
        SNVIndelReFilter_Sample(sample_flt_ann, sample_var_ann, sample_dis_ann,out_flt,out_dis1, out_dis2,json_file)

    elif tp =='Summary': 
        filelist_file=args.filelist
        flt_suffix=args.flt
        dis_suffix=args.dis
        discard_suffix=args.discard
        out_flt_suffix=args.outflt
        out_dis_suffix=args.outdis
        dirs=args.dir
        json_file=args.json
        SNVIndelReFilter_Summary(filelist_file, flt_suffix, dis_suffix, discard_suffix,out_flt_suffix, out_dis_suffix,dirs,json_file)
