# lab2mslbo

Funcionalidades para conversão dos inputs do [`SDDPlab`](https://github.com/rjmalves/sddp-lab) para pacote [`DualSDDP`](https://github.com/bfpc/DualSDDP.jl).

## Instalação

Este repositório consiste de um projeto em `julia`, não um pacote completo, de forma que inclui o arquivo `Manifest.toml` para reprodução exata do ambiente de trabalho.

Para restauração das dependências basta realizar

```
export JULIA_SSH_NO_VERIFY_HOSTS="github.com"
export JULIA_PKG_USE_CLI_GIT=true
julia
]
activate .
instantiate
```

Caso seja necessário atualizar alguma das dependências diretas de repositórios, basta fazer, no contexto do `Pkg`de `julia`, por exemplo:

```
add https://github.com/rjmalves/sddp-lab.git#branch-tag
```

## Uso

A conversão de um deck do lab se realiza através da função `build_mslbo`

```julia
include("src/lab2mslbo.jl")
modelo, dados_lab = lab2mslbo.build_mslbo("data-1dtoy/");
```

`modelo` é um objeto `DualSDDP.MSLBO` que pode ser utilizado com as demais funcionalidades do pacote e `dados_lab` é um objeto com diversas informações necessárias para compatibilização de resultados.

Mais informações estão disponíveis no arquivo `exemplo.jl`.
