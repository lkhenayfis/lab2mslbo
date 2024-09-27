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

## Uso

A conversão de um deck do lab se realiza através da função `build_mslbo`

```julia
include("src/lab2mslbo.jl")
modelo = lab2mslbo.build_mslbo("data-refactor/", seed=1234);
```

`modelo` é um objeto `DualSDDP.MSLBO` que pode ser utilizado com as demais funcionalidades do pacote.