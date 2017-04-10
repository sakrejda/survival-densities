#!/usr/local/bin/bash
expandPath() {
  local path
  local -a pathElements resultPathElements
  IFS=':' read -r -a pathElements <<<"$1"
  : "${pathElements[@]}"
  for path in "${pathElements[@]}"; do
    : "$path"
    case $path in
      "~+"/*)
        path=$PWD/${path#"~+/"}
        ;;
      "~-"/*)
        path=$OLDPWD/${path#"~-/"}
        ;;
      "~"/*)
        path=$HOME/${path#"~/"}
        ;;
      "~"*)
        username=${path%%/*}
        username=${username#"~"}
        IFS=: read _ _ _ _ _ homedir _ < <(getent passwd "$username")
        if [[ $path = */* ]]; then
          path=${homedir}/${path#*/}
        else
          path=$homedir
        fi
        ;;
    esac
    resultPathElements+=( "$path" )
  done
  local result
  printf -v result '%s:' "${resultPathElements[@]}"
  printf '%s\n' "${result%:}"
}

declare -A CONTEXT
CONTEXT_FILE=$(expandPath '~/input-store/dengue-reporting/model-weibull-uniform-mixture/processed/context.json')
export CONTEXT_FILE
CONTEXT["model-name"]=`jq -r '.["model-name"][0]' $CONTEXT_FILE`
CONTEXT["root"]=`jq -r '.root[0]' $CONTEXT_FILE`
CONTEXT["input"]=`jq -r '.input[0]' $CONTEXT_FILE`
CONTEXT["processed-input"]=`jq -r '.["processed-input"][0]' $CONTEXT_FILE`
CONTEXT["output"]=`jq -r '.output[0]' $CONTEXT_FILE`
CONTEXT["processed-output"]=`jq -r '.["processed-output"][0]' $CONTEXT_FILE`
CONTEXT["data"]=`jq -r ".data[0]" $CONTEXT_FILE`
CONTEXT["chains"]=`jq -r ".chains[0]" $CONTEXT_FILE`
export CONTEXT

