import argparse
import requests
import re


def supply_args():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('-id', help="Moon ID")
  parser.add_argument('-mode', help='info, post, or analyse')
  parser.add_argument('-snp', help='SNP VCF file')
  parser.add_argument('-cnv', help='CNV VCF file')
  parser.add_argument('-age', help='Patient Age', default="999")
  parser.add_argument('-gender', help='Patient Gender', default="Unknown")
  parser.add_argument('-consang', help='Is the patient Consanguinious?', default="false")
  parser.add_argument('-hp', help="HPO Terms", default="")
  parser.add_argument('-family', help="Family Members", default="none")
  parser.add_argument('-user', help="Email used for Moon")
  parser.add_argument('-token', help="Moon API token")
  parser.add_argument('-old_token', help="Old Moon API token")
  parser.add_argument('-server', help="Moon Server URL")
  parser.add_argument('-old_server', help="Old Moon Server for reposting")
  args = parser.parse_args()
  return args


  
def get_gender(moon_id, args):
  parameters = {"user_token": args.token, "user_email": args.user}
  
  patient_url = args.server + "/samples/" + moon_id + "/analysis.json"
  
  #print(patient_url)
  
  patient_analysis = requests.get(patient_url, params = parameters)
  
  analysis_json = patient_analysis.text
  
  #print(analysis_json)
  
  pattern = ',\"gender\":(.*?)\,'
  
  gender = re.search(pattern, analysis_json).group(1)
  
  gender_file = open("gender.txt", 'w')
  gender_file.write(gender)
  gender_file.close()
  
  return gender

def get_hpo_terms(moon_id, args):
  parameters = {"user_token": args.token, "user_email": args.user}
  
  patient_url = args.server + "/samples/" + moon_id + "/analysis.json"
  
  #print(patient_url)
  
  patient_analysis = requests.get(patient_url, params = parameters)
  
  analysis_json = patient_analysis.text
  
  #print(analysis_json)
  
  pattern = ',\"hpo_terms\":\[(.*?)\]'
  
  hpo_terms = re.search(pattern, analysis_json).group(1)
  
  hpo_file = open("hpo.txt", 'w')
  hpo_file.write(hpo_terms)
  hpo_file.close()
  
  return "hpo.txt"
  


def get_patient_info(moon_id, args):
  parameters = {"user_token": args.token, "user_email": args.user}

  patient_url = args.server + "/samples/" + moon_id + "/patient-info"
  patient_info = requests.get(patient_url, params = parameters)


  info = open("info.txt", 'w')
  

  #info.write(patient_info.content)
  info.write(patient_info.text)
  info.close()
  return "info.txt"

def post_sample(snp, args, cnv = "none", age = 999, gender = "unknown", consang = "false", hp = "", family = "none"):
  at_snp = "@/Users/campbena/Documents/galaxy-dev/tools/my_tools/" + snp
  at_cnv = "@/Users/campbena/Documents/galaxy-dev/tools/my_tools/" + cnv
  int_age = int(age)
  
  parameters = {"user_token": (None, args.token), "user_email": (None, args.user), "snp_vcf_file": (snp, open(snp, 'rb')), "sv_vcf_file": (cnv, open(cnv, 'rb')), "age": (None, age), "gender": (None, gender), "is_consanguinous": (None, consang), "hpo_terms": (None, hp), "family_members": (None, family)}
  
  if family == "none":
    parameters.pop("family_members")
    
  if cnv == "none":
    parameters.pop("sv_vcf_file")
    
  if age == "999":
    parameters.pop("age")
    
  patient_url = args.server + "/samples.json"
  
  post_to_moon = requests.post(patient_url, files = parameters)
  print(post_to_moon.text)
  print(parameters)
  
  id_file = open("new_id.txt", 'w')
  
  id_file.write(post_to_moon.text)
  id_file.close()
  
  return post_to_moon.text

def repost(args):
  #get info
  parameters = {"user_token": args.old_token, "user_email": args.user}
  
  patient_url = args.old_server + "/samples/" + args.id + "/analysis.json"
  
  patient_analysis = requests.get(patient_url, params = parameters)
  
  analysis_json = patient_analysis.text
  print(analysis_json)
  
  #get gender
  pattern = ',\"gender\":(.*?)\,'
  gender = re.search(pattern, analysis_json).group(1)
  gender = gender.replace('\"', '')
  print(gender)
  
  pattern = ',\"age\":(.*?)\,'
  age = re.search(pattern, analysis_json).group(1)
  print(age)
  
  pattern = ',\"is_consanguinous\":(.*?)\,'
  consang = re.search(pattern, analysis_json).group(1)
  print(consang)
  
  #get hpo 
  pattern = ',\"hpo_terms\":\[(.*?)\]'
  hpo_terms = re.search(pattern, analysis_json).group(1)
  hpo_terms = hpo_terms.replace('\",\"', ';')
  hpo_terms = hpo_terms.replace('\"', '')
  print(hpo_terms)
  
  return post_sample(args.snp, args, cnv = args.cnv, age = age, gender = gender, consang = consang, hp = hpo_terms, family = args.family)


def analyse_sample(moon_id, args):
  parameters = {"user_token": args.token, "user_email": args.user}
  
  patient_url = args.server + "/samples/" + moon_id + "/analysis.json"
  
  analyse_now = requests.post(patient_url, params = parameters)
  
  return analyse_now.content

def main():
  args = supply_args()
  if args.mode == "info":
    return get_patient_info(args.id, args)
  if args.mode == "post":
    return post_sample(args.snp, args, args.cnv, args.age, args.gender, args.consang, args.hp, args.family)
    
  if args.mode == "analyse":
    return analyse_sample(args.id, args)
  
  if args.mode == "hpo":
    return get_hpo_terms(args.id, args)
  
  if args.mode == "gender":
    return get_gender(args.id, args)
  
  if args.mode == "repost":
    return repost(args)

if __name__ == "__main__":
    main()