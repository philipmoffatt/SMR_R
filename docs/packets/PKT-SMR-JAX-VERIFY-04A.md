# VERIFY-04 acceptance status correction

The operator's request was to continue strict default JAX verification. Commit
`0c0152e` landed the implementation and incorrectly marked VERIFY-04 completed
after an optional native CPU exponential profile passed. That substituted a
different math implementation for the default pure-JAX acceptance path.

VERIFY-04 is reopened as active/incomplete. Its default two-day gate passes;
its default strict thirty-day gate still fails. The native CPU profile's 960
passing comparisons remain useful diagnostic evidence, not closure. VERIFY-05
must wait for the original acceptance requirement.

This follow-up changes only packet documentation and status. Existing code,
reference artifacts, numerical budgets and historical commits are preserved.
Validation: preflight passed; report-only scope and `git diff --check` passed.
