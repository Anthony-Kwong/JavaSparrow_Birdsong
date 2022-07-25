fullsong_models

bird_stats

intro_gaps <- as.data.frame( bird_stats[5])
hist(intro_gaps$avg_gap_dur)
#hist(log(intro_gaps$avg_gap_dur))

intro_dur <- as.data.frame( bird_stats[6])
hist(intro_dur$avg_song_dur)

nointro_gaps <- as.data.frame(bird_stats2[5])
hist(nointro_gaps$avg_gap_dur)

nointro_dur <- as.data.frame(bird_stats2[6])
hist(nointro_dur$avg_song_dur)
