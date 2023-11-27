import subprocess
import datetime

countries = ["slovenia", "belgium"]
weeks = ["1", "2"]

for country in countries:
    for week in weeks:
        now = datetime.datetime.now()
        print(f"{country} week {week} started at {now.hour}:{now.minute} on the {now.day}th")
        # Run BEAST
        subprocess.run(f"~/Documents/BEASTv1.10.5pre_thorney_0.1.2/bin/beast -overwrite -working gisaid_data/{country}_wk_{week}.xml", shell=True)
        # Compute annotated tree
        subprocess.run(f"~/Documents/BEASTv1.10.5pre_thorney_0.1.2/bin/treeannotator -heights ca gisaid_data/{country}_wk_{week}.trees gisaid_data/{country}_wk_{week}.nex", shell=True)