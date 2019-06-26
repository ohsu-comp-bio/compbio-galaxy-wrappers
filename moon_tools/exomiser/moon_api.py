import argparse
import requests


def supply_args():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('-id', help="Moon ID")
  parser.add_argument('-mode', help='info, post, or analyse')
  parser.add_argument('-snp', help='SNP VCF file')
  parser.add_argument('-cnv', help='CNV VCF file')
  parser.add_argument('-age', help='Patient Age')
  parser.add_argument('-gender', help='Patient Gender')
  parser.add_argument('-consang', help='Is the patient Consanguinious?')
  parser.add_argument('-hp', help="HPO Terms File")
  parser.add_argument('-family', help="Family Members File", default="none")
  args = parser.parse_args()
  return args



def get_patient_info(moon_id):
  parameters = {"user_token": "iSUQvGmVNSjq834g9fP5", "user_email": "campbena@ohsu.edu"}

  patient_url = "https://oregon.moon.diploid.com/samples/" + moon_id + "/patient-info"
  patient_info = requests.get(patient_url, params = parameters)


  info = open("info.txt", 'w')
  

  info.write(patient_info.content)
  info.close()
  return "info.txt"

def post_sample(snp, cnv, age, gender, consang, hp, family = "none" ):
  at_snp = "@/Users/campbena/Documents/galaxy-dev/tools/my_tools/" + snp
  at_cnv = "@/Users/campbena/Documents/galaxy-dev/tools/my_tools/" + cnv
  int_age = int(age)
  parameters = {"user_token": "iSUQvGmVNSjq834g9fP5", "user_email": "campbena@ohsu.edu", "snp_vcf_file": at_snp, "sv_vcf_file": at_cnv, "age": age, "gender": gender, "is_consanguinous": consang, "hpo_terms": hp, "family_members": family}
  
  
  if family == "none":
    parameters.remove("family_members")
    
  patient_url = "https://oregon.moon.diploid.com/samples.json"
  
  post_to_moon = requests.post(patient_url, params = parameters)
  print(post_to_moon.content)
  print(parameters)
  
  return post_to_moon.content

def analyse_sample(moon_id):
  parameters = {"user_token": "iSUQvGmVNSjq834g9fP5", "user_email": "campbena@ohsu.edu"}
  
  patient_url = "https://oregon.moon.diploid.com/samples/" + moon_id + "/analysis.json"
  
  analyse_now = requests.post(patient_url, params = parameters)
  
  return analyse_now.content

def main():
  args = supply_args()
  if args.mode == "info":
    return get_patient_info(args.id)
  if args.mode == "post":
    return post_sample(args.snp, args.cnv, args.age, args.gender, args.consang, args.hp, args.family)
  if args.mode == "analyse":
    return analyse_sample(args.id)

if __name__ == "__main__":
    main()
