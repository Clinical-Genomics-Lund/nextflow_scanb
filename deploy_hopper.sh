DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

LATEST_CONTAINER_BUILD="$( ls -t $DIR/container/scanb_rnaseq_*.sif |head -n1)"
CONTAINER_BASENAME=${LATEST_CONTAINER_BUILD##*/}
PIPELINE_DEST="/fs1/sima/rnaseq_scanb_dev"
CONTAINER_DEST=/fs1/resources/containers/$CONTAINER_BASENAME
DEST_HOST="rs-fs1.lunarc.lu.se"


## Deploy container if it isn't already deployed
#if test -f "$CONTAINER_DEST"; then
#    echo "Latest container already deployed, skipping!"
#else
#    echo "Deploying container"
#    scp $LATEST_CONTAINER_BUILD $DEST_HOST:$CONTAINER_DEST
#    # TODO: Replace "active" container symlink on hopper!
#fi


# Copy pipeline script
scp $DIR/main.nf $DEST_HOST:$PIPELINE_DEST

# Copy configuration file
scp $DIR/configs/nextflow.hopper.config $DEST_HOST:$PIPELINE_DEST/nextflow.config

# Copy other files
#scp -r $DIR/bin $DEST_HOST:$PIPELINE_DEST

#Added by sima:
#scp $DIR/ $DEST_HOST:$PIPELINE_DEST
#scp $DIR/run.sh $DEST_HOST:$PIPELINE_DEST
#scp ./ref/* $DEST_HOST:/fs1/resources/ref/hg38/rnaseq_scanb
git rev-parse HEAD > git.hash
scp $DIR/git.hash $DEST_HOST:$PIPELINE_DEST
