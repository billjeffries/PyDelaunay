def create_sequence_similarity_cluster_files(input_filename, pct_similar, file_limit=5):
    with open(input_filename) as f:
        clusters = f.readlines() 
        cluster_num = 1 
        for c in clusters: 
            output_filename = "data/sequence_sim_{0}_cluster_{1}.txt".format(pct_similar, cluster_num)
            with open(output_filename, "w") as cf:
                ids = c.split(' ') 
                for id in ids:
                    chain_number = int(id[5:])
                    chain_letter = chr(64+chain_number)
                    new_id = id[0:5] + chain_letter
                    cf.write(new_id + '\n')
            cluster_num += 1
            if cluster_num > file_limit:
                return

create_sequence_similarity_cluster_files("/Users/billjeffries/downloads/clusters-by-entity-30.txt", "30")

    
