shopt -s extglob
shopt -s nocaseglob

_damask_mesh() {
    local cur prev
    COMPREPLY=()

    cur=${COMP_WORDS[COMP_CWORD]}
    prev=${COMP_WORDS[COMP_CWORD-1]}
    [[ $prev == = ]] && prev=${COMP_WORDS[COMP_CWORD-2]}

    case "$prev" in
        -g|--geom|--geomfile)
            compopt -o plusdirs 2>/dev/null
            COMPREPLY=( $(compgen -G "${cur}*.msh") )
            return
            ;;
        -l|--load|--loadcase|-m|--material|--materialconfig|-n|--numerics|--numericsconfig)
            compopt -o plusdirs 2>/dev/null
            COMPREPLY=( $(compgen -G "${cur}*.@(yaml|yml)") )
            return
            ;;
        -w|--wd|--workingdir|--workingdirectory)
            COMPREPLY=( $(compgen -d -- "$cur") )
            return
            ;;
    esac

    local flags=(-g --geom --geomfile
                 -l --load --loadcase
                 -m --material --materialconfig
                 -n --numerics --numericsconfig
                 -j --job --jobname
                 -w --wd --workingdir --workingdirectory
                 -h --help)
    if (( COMP_CWORD == 1 )) || [[ $cur == -* ]] || [[ -z $cur ]]; then
        COMPREPLY=($(compgen -W "${flags[*]}" -- "$cur"))
        compopt -o nosort 2>/dev/null
        return
    fi
}

complete -F _damask_mesh DAMASK_mesh
