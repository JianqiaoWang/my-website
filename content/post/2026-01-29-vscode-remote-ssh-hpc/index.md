---
title: "VS Code、Remote-SSH、高校集群与计算环境"
date: 2026-01-29
author: Jianqiao Wang
tags:
  - HPC
  - Remote-SSH
  - VSCode
  - 科研计算
categories:
  - 技术随笔
---

## VS Code、Remote-SSH、集群以及计算环境

最近在身边一些“技术先锋”的带动下，我逐渐意识到自己的科研工具链似乎有些跟不上时代。于是，这段时间我开始主动学习和尝试新工具，探索 AI coding，希望找到一种能够将其合理融入现有科研流程的平衡点。在深度体验了 Cursor、Positron 以及各类 AI 插件之后，我逐渐形成了一个较为强烈的主观感受：目前的 AI coding 更多是工程驱动的，而非为数据科学或科研问题本身而设计。它更像是一个高效的“需求实现者”——只要问题被清晰地定义，它往往能给出不错的实现方案；但它并不擅长帮助研究者探索和重塑问题本身。而科研恰恰是一个不断试探、反复迭代、持续修正研究问题的过程，最终形成的问题往往与最初的直觉相去甚远。在这种探索过程中，AI 对代码的频繁修改反而容易破坏原本简洁清晰的结构，面对被“改乱”的代码，难免让人感到有些沮丧。

当然，如何让 AI 更好地服务科研工作者，是一个非常宏大且尚未有定论的问题。本文并不试图展开这一讨论，而是聚焦于一个更具体、更工具层面的话题：新工具的快速演进，与高校计算基础设施之间逐渐显现出的错位。

在我的大多数科研工作中，核心数据与程序并不存放在本地，而是托管在学校的计算集群上，因此本地与远端之间的连接始终是刚需。过去，我相对稳定的工作流主要依赖“老三样”：**MobaXterm（上传文件、提交 sbatch 作业） + 本地 RStudio（代码编辑） + Web 界面（VDI 调试、查资料或使用 ChatGPT）**。

但另一方面，诸如 VS Code、Cursor、Positron 等工具已经高度集成了插件生态和 AI 辅助能力。例如 Remote-SSH 可以显著提升远程开发的效率，而 Codex、CC、Copilot 等插件能够直接理解代码库结构，针对性地生成或补全代码。如果这些工具能够顺利接入集群环境，许多重复而繁琐的工作本应可以被大幅简化。

因此，在过去一段时间里，我系统性地尝试使用 VS Code 的 Remote-SSH 功能连接不同高校的高性能计算集群（HPC），包括清华的探索1000、UPenn 的 PMACS，以及哈佛的 FASRC。过程中踩了不少坑，也逐渐形成了一些相对稳定的认识。这里做一个阶段性的整理，一方面作为个人记录，另一方面也希望能为有类似需求的科研用户提供一些参考。 

### Remote-SSH 在高校集群中的结构性限制

最初在配置 Remote-SSH 时，我一度以为问题出在 VS Code 或具体参数设置上。但在反复尝试后逐渐意识到，在高校集群环境中，Remote-SSH 的困难往往并非工具层面的，而是由集群的整体结构所决定。在大多数集群中，用户只能直接登录 **login node**，而真正承担计算与开发负载的 **compute node** 则由调度系统动态分配、默认不可直连。这一问题本身，已经属于集群设计层面的取舍，而非用户侧可以通过技巧解决。

在实际使用中，真正让我反复受阻的，并不是调度系统本身，而是操作系统层面的限制。清华的探索1000集群以及 PMACS HPC 至今仍运行在 **CentOS 7.9**。这一选择在 2025 年已经很难再为其合理性辩护。随着开发工具的快速迭代，老旧系统逐渐成为瓶颈。例如，新版本的 VS Code 已不再支持 CentOS 7.9，这意味着在探索1000集群中，用户只能被迫长期停留在非常旧的 VS Code 版本，几乎没有升级空间。

在对比之下，哈佛的 FASRC 显得尤为成熟。他们已整体升级至 Rocky Linux 8，这一决策本身表明，管理者清楚地认识到科研开发环境是不断演进的，系统应当为新工具预留空间，而不是要求用户长期迁就旧平台。与此同时，FASRC 也认真对待“从本地进行远程开发”的实际需求，明确支持通过 SSH、Remote-SSH 或 Tunnel 等方式连接计算节点，并将这些路径写入官方文档进行维护。相关说明可以在其官方文档中找到： [https://docs.rc.fas.harvard.edu/kb/vscode-remote-development-via-ssh-or-tunnel/](https://docs.rc.fas.harvard.edu/kb/vscode-remote-development-via-ssh-or-tunnel/)

当然受限于集群这一个大框架，他们并不建议用vscode直连。相反

> We encourage all our users to utilize the Open On Demand (OOD) web interface of the cluster to launch VS Code when remote development work is not required. The instructions to launch VS Code using the Remote Desktop App are given [here](https://docs.rc.fas.harvard.edu/kb/ood-remote-desktop-how-to-open-software/#Visual_Studio_Code).

这段官方说明也很有代表性：在无法或不需要进行真正的远程开发时，FASRC 明确引导用户使用基于 Web 的 OOD/远程桌面方案，而不是强行将所有工作流压到 SSH 之上。

PMACS HPC 同样提供了较为成熟的远程桌面方案，而清华探索1000目前尚未提供类似支持。（THU计算机底子还是厚，不需要这种傻瓜式操作）

当在经历了一系列尝试之后，我最终找到了一条相对可行、且个人实践中效果还挺不错的解决路径，其具体技术方案在下面的仓库中，供有需要的读者参考：

[https://github.com/MikeWang000000/vscode-server-centos7?tab=readme-ov-file](https://github.com/MikeWang000000/vscode-server-centos7?tab=readme-ov-file)

最后，我仍然想说许多现有的计算环境更多是为传统工程或工科场景设计的，对以数据分析为核心的新型科研范式支持仍然有限。在数据科学领域，科研工作已高度依赖 IDE 驱动的交互式开发过程——研究者需要在代码、数据与结果之间频繁切换，不断调试、可视化并重构分析流程。研究过程中，判断分析路径是否合理、流程是否简洁、结果是否可信，往往比单次计算本身更为关键。实际上，我常常前进一步，就要后退两步检查一下。如果一个系统仍然默认科研计算只是“写脚本—交作业—等结果”的线性流程，那么它在客观上就难以为这种现代研究方式提供足够的支持。也希望单位的相关基础设施能够逐步向这一方向靠拢。
