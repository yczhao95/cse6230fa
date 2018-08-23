#!/bin/bash

if [ -z $GTUSER ]; then
  echo "Define GTUSER to be your GT user name (e.g. gburdell3) to use this script"
  exit 1
fi

if [ -z $1 ]; then
  echo "Usage: $0 BRANCH ASSIGNMENT_NAME"
  echo "Give the name of a branch or a commit:"
  echo "this will become the master branch of the assignment repository"
  exit 1
fi

if [ -z $2 ]; then
  echo "Usage: $0 BRANCH ASSIGNMENT_NAME"
  echo "Give the name of the assignment (e.g. as1):"
  echo "this will become the master branch of the assignment repository"
  exit 1
fi

REV_NAME="$1"
REV_PARSE=`git rev-parse --short $1` || exit $?
ASSIGN_NAME="$2"

REMOTE_NAME="cse6230-$ASSIGN_NAME-$GTUSER"

ASSIGN_REMOTE="https://$GTUSER@github.gatech.edu/$GTUSER/$REMOTE_NAME.git"

echo "Pushing revision $REV_NAME ($REV_PARSE) to repository https://github.gatech.edu/$GTUSER/${REMOTE_NAME}.git"

# https://stackoverflow.com/questions/1885525/how-do-i-prompt-a-user-for-confirmation-in-bash-script
read -p "Continue? (y/n)" -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi

if `git remote | grep --quiet ${ASSIGN_NAME}` ; then
  read -p "Remote ${ASSIGN_NAME} already exists, are you sure? (y/n)" -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]
  then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
  fi
else
  echo "Creating remote ${ASSIGN_NAME} (${ASSIGN_REMOTE})"
  git remote add $ASSIGN_NAME $ASSIGN_REMOTE || exit $?
  HTTP_CODE=`curl -s -u $GTUSER https://github.gatech.edu/api/v3/user/repos -d "{\"name\":\"$REMOTE_NAME\",\"private\":true}" -o /dev/null -w "%{http_code}"`
  if [ ! "x$HTTP_CODE" = "x201" ]; then
    echo "Failed to creat repo: return code $HTTP_CODE"
    exit 1
  fi
fi

echo "Pushing $REV_NAME ($REV_PARSE) to ${ASSIGN_NAME}"
git push $ASSIGN_NAME $REV_PARSE:refs/heads/master 
echo "Transfering the repository to cse6230fa18"
HTTP_CODE=`curl -s -u $GTUSER https://github.gatech.edu/api/v3/repos/${GTUSER}/${REMOTE_NAME}/transfer -H "Accept: application/vnd.github.nightshade-preview+json" -d "{\"new_owner\":\"cse6230fa18\",\"team_ids\":[3994]}" -o /dev/null -w "%{http_code}"`
if [ ! "x$HTTP_CODE" = "x202" ]; then
  echo "Failed to creat repo: return code $HTTP_CODE"
  exit 1
fi
echo "Transfer complete!"
