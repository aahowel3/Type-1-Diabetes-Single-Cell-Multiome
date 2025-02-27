## Downloading TopMed imputation results
# Unzip and check files
# Replace qvMG1IL@Gq0gj with the password in the topmed email
for file in *.zip; do unzip -P qvMG1IL@Gq0gj $file; done
md5sum *.zip
