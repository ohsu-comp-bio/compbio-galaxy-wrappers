PROJECT_ID=$1

eval $(bs load config ilssd)
bs appresult download -i ${PROJECT_ID}
