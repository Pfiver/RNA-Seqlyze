
git checkout 4efe4317c377019fbb64bcfae519403a3a9d1f5f

make -C bcftools && make -C misc && make SUBDIRS=. LIBCURSES= DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE" "AOBJS=bam_plcmd.o sam_view.o bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o cut_target.o phase.o bam2depth.o"

for d in bam lib include; do ln -s . $d; done
