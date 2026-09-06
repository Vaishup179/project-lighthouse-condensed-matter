# Project Lighthouse Validation

Local validation performed on 6 September 2026 using Harbor with Docker Desktop through Ubuntu WSL2.

## Oracle validation

Command:

```bash
harbor run -p tasks/physics/condensed-matter-physics/sdh-pocket-reconstruction -a oracle
```

Terminal output:

```text
  1/1 Mean: 1.000 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 0:00:43 0:00:00
adhoc • oracle
┏━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━┓
┃ Trials ┃ Exceptions ┃  Mean ┃
┡━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━┩
│      1 │          0 │ 1.000 │
└────────┴────────────┴───────┘

┏━━━━━━━━┳━━━━━━━┓
┃ Reward ┃ Count ┃
┡━━━━━━━━╇━━━━━━━┩
│ 1.0    │     1 │
└────────┴───────┘

Job Info
Total runtime: 45s
Results written to jobs/2026-09-06__15-51-31/result.json
Inspect results by running `harbor view jobs`
Share results by running `harbor upload jobs/2026-09-06__15-51-31`

```

**Oracle mean: 1.000**

## NOP validation

Command:

```bash
harbor run -p tasks/physics/condensed-matter-physics/sdh-pocket-reconstruction -a nop
```

Terminal output:

```text
  1/1 Mean: 0.000 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 0:00:10 0:00:00
adhoc • nop
┏━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━┓
┃ Trials ┃ Exceptions ┃  Mean ┃
┡━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━┩
│      1 │          0 │ 0.000 │
└────────┴────────────┴───────┘

┏━━━━━━━━┳━━━━━━━┓
┃ Reward ┃ Count ┃
┡━━━━━━━━╇━━━━━━━┩
│ 0.0    │     1 │
└────────┴───────┘

Job Info
Total runtime: 11s
Results written to jobs/2026-09-06__15-53-23/result.json
Inspect results by running `harbor view jobs`
Share results by running `harbor upload jobs/2026-09-06__15-53-23`

```

**NOP mean: 0.000**
