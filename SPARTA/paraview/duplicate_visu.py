h_s = [100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]

content_ref = open("visu_%ikm.pvsm" % h_s[0]).read()
for h in h_s[1:]:
    content_new = content_ref.replace("MAV_stage_2_%ikm"%h_s[0], "MAV_stage_2_%ikm"%h)
    content_new = content_new.replace("h=%ikm"%h_s[0], "h=%ikm"%h)
    f = open("visu_%ikm.pvsm" % h, "w")
    f.write(content_new)
    f.close()